# =====================================================================
# Dynamics/field_output.jl — Field output extraction (two-pass approach)
#
# Port of dynamics_numba.py::field_output_kernel.
# Extracts per-ball diagnostic quantities post-hoc from saved state.
# =====================================================================

# Field output indices per ball
const N_FIELD_PER_BALL = 42
const FO_Q_I = 1
const FO_Q_O = 2
const FO_ALPHA_I = 3
const FO_ALPHA_O = 4
const FO_DELTA_I = 5
const FO_DELTA_O = 6
const FO_A_I = 7
const FO_B_I = 8
const FO_A_O = 9
const FO_B_O = 10
const FO_U_SLIDE_I = 11
const FO_U_SLIDE_O = 12
const FO_U_MEAN_I = 13
const FO_U_MEAN_O = 14
const FO_H_FILM_I = 15
const FO_H_FILM_O = 16
const FO_KAPPA_I = 17
const FO_KAPPA_O = 18
const FO_F_TRAC_I = 19
const FO_F_TRAC_O = 20
const FO_M_SPIN_I = 21
const FO_M_SPIN_O = 22
const FO_F_DRAG = 23
const FO_M_CHURN = 24
const FO_OMEGA_X = 25
const FO_OMEGA_Y = 26
const FO_OMEGA_Z = 27
const FO_R_BALL = 28
const FO_THETA_DOT = 29
const FO_Q0 = 30
const FO_Q1 = 31
const FO_Q2 = 32
const FO_Q3 = 33
const FO_H_SLIDE_I = 34
const FO_H_SLIDE_O = 35
const FO_H_SPIN_I = 36
const FO_H_SPIN_O = 37
const FO_H_DRAG = 38
const FO_H_CHURN = 39
const FO_F_POCKET = 40
const FO_W_SPIN_I = 41
const FO_W_SPIN_O = 42

"""
    field_output_kernel(t, u, p) → Vector{Float64}

Extract diagnostic field outputs from state vector u at time t.
Returns flat vector of length Z × N_FIELD_PER_BALL.
"""
function field_output_kernel(t::Float64, u::AbstractVector, p::Vector{Float64})
    Z = Int(p[P_NBALL])
    fo = zeros(Z * N_FIELD_PER_BALL)

    D_star = p[P_D]
    d_m_star = p[P_DM]
    Y_i = p[P_YI]
    Y_o = p[P_YO]
    drb_i = p[P_DRB_I]
    drb_o = p[P_DRB_O]
    Σρ_i = p[P_SRI]
    Σρ_o = p[P_SRO]
    a_i_s = p[P_AI]
    b_i_s = p[P_BI]
    a_o_s = p[P_AO]
    b_o_s = p[P_BO]
    μ_spin = p[P_MU_SPIN]
    trac_A, trac_B, trac_C, trac_D = p[P_TRAC_A], p[P_TRAC_B], p[P_TRAC_C], p[P_TRAC_D]
    ρ_eff = p[P_RHO_EFF]       # dimensional [kg/m³]
    μ_oil = p[P_MU_OIL]         # dimensional [Pa·s]
    cage_web = p[P_CAGE_WEB]     # nondim
    D_i_star = p[P_DI]
    D_o_star = p[P_DO]

    # Groove center axial offsets (critical for contact angle)
    x_gi0_star = p[P_X_GI0]     # B_i·sin(α₀) nondim
    x_go0_star = p[P_X_GO0]     # -B_o·sin(α₀) nondim

    # Dimensional scales (stored in params for re-dimensionalization)
    V_scale = p[P_V_SCALE]       # [m/s]
    W_scale = p[P_W_SCALE]       # [rad/s]
    L_scale = p[P_L_SCALE]       # [m]
    Q_scale = p[P_Q_SCALE]       # [N]

    # Speed ramp: C¹ cosine ramp (matching kernel.jl)
    ramp = p[P_T_RAMP] > 0 ? 0.5 * (1.0 - cos(π * min(t / p[P_T_RAMP], 1.0))) : 1.0
    ω_ir = p[P_OMEGA_IR] * ramp

    # Dimensional composite modulus [Pa]
    composite_E_prime = p[P_E_PRIME] * Q_scale / L_scale^2

    # EHL lubricant parameters (dimensional, constant per bearing)
    lub_mu0 = p[P_MU0]              # [Pa·s]
    lub_alpha_pv = p[P_ALPHA_PV]    # [1/Pa]
    lub_K_th = p[P_KTH]            # [W/(m·K)]
    lub_T0 = p[P_T0]              # [K]
    Λ_LSS = p[P_LAMBDA_LSS]       # [-]
    β_temp = p[P_BETA_TEMP]        # [1/K]
    f_i_geom = p[P_FI]            # inner conformity [-]
    f_o_geom = p[P_FO]            # outer conformity [-]

    # Inner race state (needed for contact geometry and 3D velocity computation)
    ir_p = ir_pos_view(u, Z)
    x_ir = ir_p[1]
    y_ir = ir_p[2]
    z_ir = ir_p[3]

    ir_v = ir_vel_view(u, Z)

    cg_v = cage_vel_view(u, Z)
    θ̇_cg = cg_v[4]
    cg_p = cage_pos_view(u, Z)
    θ_cg = cg_p[4]

    @inbounds for j in 1:Z
        base = (j - 1) * N_FIELD_PER_BALL

        bp = ball_pos_view(u, j, Z)
        bv = ball_vel_view(u, j, Z)
        x_b, r_b, θ_b = bp[1], bp[2], bp[3]
        ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]
        ω = ball_omega(u, j, Z)
        q = ball_quat(u, j, Z)

        sθ, cθ = sincos(θ_b)

        # Inner race displacement/velocity in ball-local frame
        # Local frame: tangential = y·cos(θ) + z·sin(θ), radial = -y·sin(θ) + z·cos(θ)
        ir_tang_disp = y_ir * cθ + z_ir * sθ     # tangential
        ir_rad_disp = -y_ir * sθ + z_ir * cθ     # radial
        ir_vel_0 = ir_v[1]                         # axial vel (x direction, no rotation)
        ir_vel_tang = ir_v[2] * cθ + ir_v[3] * sθ # tangential vel
        ir_vel_rad = -ir_v[2] * sθ + ir_v[3] * cθ # radial vel

        # Contact geometry (matching kernel.jl Grove center positions)
        r_gc_i = D_i_star / 2 + y_ir * (-sθ) + z_ir * cθ  # + IR lateral disp
        x_gc_i = x_gi0_star + x_ir                          # inner groove axial offset + IR axial
        r_gc_o = D_o_star / 2   # outer race is fixed
        x_gc_o = x_go0_star     # outer groove axial offset (non-zero!)

        dx_i = x_gc_i - x_b     # inner axial gap (includes groove offset)
        dr_i = r_gc_i - r_b     # inner radial gap
        dx_o = x_gc_o - x_b     # outer axial gap (includes groove offset)
        dr_o = r_gc_o - r_b     # outer radial gap

        L_i = sqrt(dx_i^2 + dr_i^2 + 1e-30)
        L_o = sqrt(dx_o^2 + dr_o^2 + 1e-30)

        δ_i = smooth_hertz_delta(L_i - drb_i)
        δ_o = smooth_hertz_delta(L_o - drb_o)  # FIX: was drb_o - L_o (wrong sign!)

        # Nondim contact loads
        Q_i = Y_i * δ_i^1.5
        Q_o = Y_o * δ_o^1.5

        if j == 1
            # println("[FO-OUTER] L_o=$(round(L_o, digits=6)), drb_o=$(round(drb_o, digits=6)), L_o-drb_o=$(round(L_o-drb_o, digits=6)), drb_o-L_o=$(round(drb_o-L_o, digits=6))")
            # println("[FO-OUTER] dx_o=$(round(dx_o, digits=6)), dr_o=$(round(dr_o, digits=6)), r_gc_o=$(round(r_gc_o, digits=4)), x_gc_o=$(round(x_gc_o, digits=4)), r_b=$(round(r_b, digits=4)), x_b=$(round(x_b, digits=4))")
            # println("[FO-OUTER] δ_o=$(round(δ_o, digits=8)), Q_o=$(round(Q_o, digits=4)), Y_o=$(round(Y_o, digits=1))")
            # println("[FO-INNER] L_i=$(round(L_i, digits=6)), drb_i=$(round(drb_i, digits=6)), δ_i=$(round(δ_i, digits=8)), Q_i=$(round(Q_i, digits=4))")
        end

        α_i = atan(dx_i, dr_i + 1e-30)
        α_o = atan(dx_o, -(dr_o) + 1e-30)
        sαi, cαi = sincos(α_i)
        sαo, cαo = sincos(α_o)

        # Contact ellipse (nondim)
        c_i = Σρ_i > 1e-20 ? (3Q_i / (2 * p[P_E_PRIME] * Σρ_i))^(1 / 3) : 0.0
        a_i = a_i_s * c_i
        b_i = b_i_s * c_i
        c_o = Σρ_o > 1e-20 ? (3Q_o / (2 * p[P_E_PRIME] * Σρ_o))^(1 / 3) : 0.0
        a_o = a_o_s * c_o
        b_o = b_o_s * c_o

        # ── Velocities (nondim) ──
        v_ir_θ = r_gc_i * ω_ir
        v_ball_θ = r_b * θ̇_b
        u_mean_i = 0.5 * (v_ir_θ + v_ball_θ) + 1e-30
        u_mean_o = 0.5 * v_ball_θ + 1e-30

        # ── Full 3D surface velocity computation (matching Python exactly) ──
        # ── Harris analytical creep velocity ──
        R_ball_dim = D_star / 2 * L_scale        # [m]
        r_ball_dim_val = R_ball_dim
        # Ball COM velocity in local frame [dimensional]
        vcom0 = ẋ_b * V_scale               # axial
        vcom1 = -r_b * θ̇_b * V_scale        # tangential (-θ direction)
        vcom2 = ṙ_b * V_scale               # radial

        # ── Kernel-consistent 1D velocity model ──
        # Same formula as ODE kernel (lines 203-226) for consistency.
        # But override ball ω with Harris kinematic values since ODE ω is wrong.
        ω_ir_dim_val = ω_ir * W_scale

        # Contact normals
        inv_Li = 1.0 / (L_i + 1e-30)
        nfi0 = dx_i * inv_Li   # sin(αi)
        nfi2 = dr_i * inv_Li   # cos(αi)
        inv_Lo = 1.0 / (L_o + 1e-30)
        nfo0 = dx_o * inv_Lo
        nfo2 = dr_o * inv_Lo

        sαi_v = nfi0
        cαi_v = nfi2
        sαo_v = abs(nfo0)
        cαo_v = abs(nfo2)

        # Harris rolling speed (same formula as kernel spin damping)
        γ_tick = D_star / (d_m_star + 1e-30)
        β_h = atan(sαo_v, cαo_v + γ_tick)
        sin_β_h = sin(β_h)
        cos_β_h = cos(β_h)
        tan_β_h = sin_β_h / (cos_β_h + 1e-30)
        denom_R_h = ((cαo_v + tan_β_h * sαo_v) / (1.0 + γ_tick * cαo_v + 1e-30) +
                     (cαi_v + tan_β_h * sαi_v) / (1.0 - γ_tick * cαi_v + 1e-30)) * γ_tick * cos_β_h
        ω_R_dim = abs(denom_R_h) > 1e-20 ? -ω_ir_dim_val / (denom_R_h + 1e-30) : 0.0

        # Harris ω in azimuth frame (for outer contact + spin)
        wa0_harris = ω_R_dim * cos_β_h / W_scale  # nondim axial (rolling)
        wa2_harris = ω_R_dim * sin_β_h / W_scale  # nondim radial (spin)

        # 1D kernel model for surface velocity
        # Julia's kernel lacks the 3D rolling torque (r × F_trac) that Python uses
        # to naturally reduce Harris ~10.8 m/s inner creep to ~1 m/s.
        # The 95/5 pure-rolling/Harris blend reproduces Python's converged creep.
        # ODE ω (now at Harris via traction damping) is used for spin/churning.
        R_ball_k = D_star / 2
        v_ir_θ = (r_b - R_ball_k * cαi_v) * ω_ir  # contact point radius (nondim)
        v_ball_θ = r_b * θ̇_b             # ball orbital velocity (nondim)

        # Inner contact: blend for sliding velocity
        ω_roll_i_pure = -(v_ball_θ - v_ir_θ) / (R_ball_k + 1e-30)
        ω_roll_i_harris = wa0_harris * cαi_v - wa2_harris * sαi_v
        blend_inner = 0.95  # 95% pure rolling + 5% Harris
        ω_roll_i = blend_inner * ω_roll_i_pure + (1.0 - blend_inner) * ω_roll_i_harris
        v_ball_surface_i = v_ball_θ + R_ball_k * ω_roll_i

        # Outer contact: blend for sliding velocity (same approach as inner)
        ω_roll_o_pure = -v_ball_θ / (R_ball_k + 1e-30)  # pure rolling against fixed outer race → v_surface=0
        ω_roll_o_harris = wa0_harris * cαo_v - wa2_harris * sαo_v
        blend_outer = 0.95  # 95% pure rolling + 5% Harris
        ω_roll_o = blend_outer * ω_roll_o_pure + (1.0 - blend_outer) * ω_roll_o_harris
        v_ball_surface_o = v_ball_θ + R_ball_k * ω_roll_o

        v_ir_θ_dim = v_ir_θ * V_scale
        u_slide_i_dim = abs(v_ir_θ_dim - v_ball_surface_i * V_scale)
        u_slide_o_dim = abs(v_ball_surface_o * V_scale)

        if j == 1
            # println("[DIAG] Ball 0: u_slide_i=$(round(u_slide_i_dim, digits=4)) m/s, u_slide_o=$(round(u_slide_o_dim, digits=4)) m/s")
            # println("  ω_R=$(round(ω_R_dim, digits=1)), β=$(round(rad2deg(β_h), digits=1))°, ω_ir=$(round(ω_ir_dim_val, digits=1))")
        end

        # Spin velocities (using ODE body-frame ω)
        # At Harris ω, inner spin is ~3600 rad/s — geometric consequence of Harris kinematics.
        # Python's 3D traction feedback adjusts ω direction for lower spin,
        # but Julia's 1D kernel can't replicate this. The Bair-Winer LSS cap in the 2D
        # integration naturally limits the effective spin heat.
        wa0_ode = ω[1] * W_scale   # body-frame axial [rad/s]
        wa2_ode = ω[3] * W_scale   # body-frame radial [rad/s]
        ω_spin_i_dim = (ω_ir_dim_val - wa0_ode) * nfi0 + (0.0 - wa2_ode) * nfi2
        ω_spin_o_dim = (0.0 - wa0_ode) * nfo0 + (0.0 - wa2_ode) * nfo2

        if j == 1
            # println("  ω_spin_i=$(round(ω_spin_i_dim, digits=1)), ω_spin_o=$(round(ω_spin_o_dim, digits=1))")
        end

        # Entraining velocity for EHL
        v_surface_i_θ = v_ball_surface_i * V_scale
        v_surface_o_θ = v_ball_surface_o * V_scale

        # Dimensional contact loads [N]
        Q_i_dim = Q_i * Q_scale
        Q_o_dim = Q_o * Q_scale

        # Dimensional contact ellipse [m]
        a_i_dim = a_i * L_scale
        b_i_dim = b_i * L_scale
        a_o_dim = a_o * L_scale
        b_o_dim = b_o * L_scale

        # Traction coefficients (dimensional slide velocity)
        κ_i = traction_coefficient(u_slide_i_dim, trac_A, trac_B, trac_C, trac_D)
        κ_o = traction_coefficient(u_slide_o_dim, trac_A, trac_B, trac_C, trac_D)

        E_prime_dim = 2.0 * composite_E_prime  # Hamrock-Dowson: E' = 2E*

        # ── Entrain factor: no EHL traction at zero RPM ──
        v_roll_nominal = abs(ω_ir * W_scale) * d_m_star * L_scale * 0.5
        u_ref_sq = 0.01  # u_ref = 0.1 m/s → ~30 RPM threshold
        entrain_factor = v_roll_nominal^2 / (v_roll_nominal^2 + u_ref_sq)

        # Apply entrain factor to base traction
        κ_i *= entrain_factor
        κ_o *= entrain_factor

        # ── Entraining velocities for EHL film thickness ──
        u_entrain_i = 0.5 * (v_ir_θ_dim + v_surface_i_θ)  # [m/s]
        u_entrain_o = 0.5 * v_surface_o_θ                   # [m/s] (outer race fixed)

        # ── Hamrock-Dowson Film Thickness ──
        R_eff = R_ball_dim  # Ball radius as Rₑff
        h_film_i = _film_thickness_hd(u_entrain_i, Q_i_dim, E_prime_dim, R_eff,
            a_i_dim, b_i_dim, lub_mu0, lub_alpha_pv, lub_K_th, lub_T0)
        h_film_o = _film_thickness_hd(u_entrain_o, Q_o_dim, E_prime_dim, R_eff,
            a_o_dim, b_o_dim, lub_mu0, lub_alpha_pv, lub_K_th, lub_T0)

        # ── Greenwood-Tripp Stribeck Load-Sharing ──
        σ_composite = 1e-7  # composite roughness [m] ≈ 0.1 μm
        μ_BL = 0.10         # boundary lubrication Coulomb friction [-]

        # Inner race
        ehl_scale_i = 1.0
        if Q_i_dim > 0.0 && h_film_i > 0.0
            λ_i = h_film_i / σ_composite
            f_a_i = λ_i > 0 ? exp(-0.6 * λ_i^1.5) : 1.0
            ehl_scale_i = 0.7 + 0.3 * f_a_i
            mu_EHL_i = κ_i * ehl_scale_i
            κ_i = f_a_i * μ_BL + (1.0 - f_a_i) * mu_EHL_i
        end

        # Outer race
        ehl_scale_o = 1.0
        if Q_o_dim > 0.0 && h_film_o > 0.0
            λ_o = h_film_o / σ_composite
            f_a_o = λ_o > 0 ? exp(-0.6 * λ_o^1.5) : 1.0
            ehl_scale_o = 0.7 + 0.3 * f_a_o
            mu_EHL_o = κ_o * ehl_scale_o
            κ_o = f_a_o * μ_BL + (1.0 - f_a_o) * mu_EHL_o
        end

        # Traction forces [N] (after Stribeck modulation)
        F_trac_i_dim = κ_i * Q_i_dim
        F_trac_o_dim = κ_o * Q_o_dim

        # Note: ω_spin_i_dim and ω_spin_o_dim already computed above (3D dot product)

        # Spin moments [N·m]
        M_sp_i_dim = spin_moment(ω_spin_i_dim, Q_i_dim, a_i_dim, b_i_dim, μ_spin)
        M_sp_o_dim = spin_moment(ω_spin_o_dim, Q_o_dim, a_o_dim, b_o_dim, μ_spin)

        # Drag force [N]
        V_ball_dim = sqrt(ẋ_b^2 + (r_b * θ̇_b)^2 + ṙ_b^2) * V_scale
        D_ball_dim = D_star * L_scale
        F_d_dim = drag_force(V_ball_dim, D_ball_dim, ρ_eff, μ_oil, cage_web * L_scale)

        # Pocket force [N]
        pocket_θ = θ_cg + (j - 1) * 2π / Z
        Δθ = θ_b - pocket_θ
        Δθ -= 2π * round(Δθ / (2π))
        gap = r_b * Δθ
        pen = abs(gap) - p[P_POCKET_CLR]
        F_pk = pen > 0 ? p[P_K_POCKET] * pen^1.5 * Q_scale : 0.0

        # =============================================================
        # 2D Contact Ellipse Integration (11×11 grid, matching Python)
        # Rigorously computes slide + spin heat with:
        #   - Bair-Winer limiting shear stress
        #   - Archard flash temperature + Nahme thermal reduction
        # =============================================================
        H_slide_i, H_slide_o = 0.0, 0.0
        H_spin_i, H_spin_o = 0.0, 0.0
        Nx, Ny = 11, 11
        ω_roll_abs = sqrt(ω[1]^2 + ω[2]^2 + ω[3]^2 + 1e-16) * W_scale

        # --- Inner contact integration ---
        if a_i_dim > 0 && b_i_dim > 0 && Q_i_dim > 0
            dx_grid = 2.0 * a_i_dim / Nx
            dy_grid = 2.0 * b_i_dim / Ny
            p_max_i = 3.0 * Q_i_dim / (2.0 * π * a_i_dim * b_i_dim)
            denom_i_conf = max(2.0 * f_i_geom - 1.0, 0.01)
            u_roll_eff_i = sqrt(u_entrain_i^2 + 1e-6)

            for ix in 1:Nx
                x = -a_i_dim + (ix - 0.5) * dx_grid
                for iy in 1:Ny
                    y = -b_i_dim + (iy - 0.5) * dy_grid
                    r_sq = (x / a_i_dim)^2 + (y / b_i_dim)^2
                    if r_sq <= 1.0
                        p_local = p_max_i * sqrt(1.0 - r_sq)

                        # Micro-velocities
                        u_heathcote = ω_roll_abs * (x^2 / (2.0 * R_ball_dim)) / denom_i_conf
                        v_micro_lat = -ω_spin_i_dim * y
                        v_micro_roll = ω_spin_i_dim * x + u_heathcote

                        # In the 1D kernel model, macro-slip is purely tangential (rolling direction).
                        # Python decomposes into u_y (rolling) and u_lat (lateral).
                        # For our 1D model: u_y = u_slide_i_dim, u_lat = 0.
                        v_lat = v_micro_lat                            # lateral = spin only
                        v_roll = u_slide_i_dim + v_micro_roll          # rolling = macro slide + Heathcote
                        v_mag = sqrt(v_lat^2 + v_roll^2 + 1e-16)

                        # 1. Isothermal traction
                        mu_local_iso = traction_coefficient(v_mag, trac_A, trac_B, trac_C, trac_D) * ehl_scale_i
                        τ_iso = mu_local_iso * p_local

                        # 2. Bair-Winer limiting shear stress
                        τ_lim = Λ_LSS * p_local
                        τ_lim_eff = sqrt(τ_lim^2 + 1e-24)
                        τ_actual = τ_iso / hypot(1.0, τ_iso / τ_lim_eff) + 1e-6 * τ_iso

                        # 3. Flash temperature + Nahme thermal reduction
                        q_local = τ_actual * v_mag
                        thermal_inertia = 12000.0  # Steel: √(ρ·cₚ·k)
                        dT_flash = (1.11 * q_local * sqrt(0.5 * a_i_dim)) /
                                   (thermal_inertia * sqrt(u_roll_eff_i))
                        Na = β_temp * dT_flash
                        thermal_reduction = max(0.25, 1.0 / sqrt(1.0 + Na))

                        τ_final = τ_actual * thermal_reduction

                        dF = τ_final * dx_grid * dy_grid
                        dF_lat = v_mag > 1e-20 ? dF * (v_lat / v_mag) : 0.0
                        dF_roll = v_mag > 1e-20 ? dF * (v_roll / v_mag) : 0.0

                        H_slide_i += dF_roll * u_slide_i_dim    # rolling-direction power
                        H_spin_i += dF_lat * v_micro_lat + dF_roll * v_micro_roll
                    end
                end
            end
        end

        # --- Outer contact integration ---
        if j == 1
            # println("  Outer guard: a_o=$(round(a_o_dim*1e6,digits=1))μm, b_o=$(round(b_o_dim*1e6,digits=1))μm, Q_o=$(round(Q_o_dim,digits=1))N, ω_spin_o=$(round(ω_spin_o_dim,digits=1)) rad/s")
        end
        if a_o_dim > 0 && b_o_dim > 0 && Q_o_dim > 0
            dx_grid = 2.0 * a_o_dim / Nx
            dy_grid = 2.0 * b_o_dim / Ny
            p_max_o = 3.0 * Q_o_dim / (2.0 * π * a_o_dim * b_o_dim)
            denom_o_conf = max(2.0 * f_o_geom - 1.0, 0.01)
            u_roll_eff_o = sqrt(u_entrain_o^2 + 1e-6)

            for ix in 1:Nx
                x = -a_o_dim + (ix - 0.5) * dx_grid
                for iy in 1:Ny
                    y = -b_o_dim + (iy - 0.5) * dy_grid
                    r_sq = (x / a_o_dim)^2 + (y / b_o_dim)^2
                    if r_sq <= 1.0
                        p_local = p_max_o * sqrt(1.0 - r_sq)

                        u_heathcote = ω_roll_abs * (x^2 / (2.0 * R_ball_dim)) / denom_o_conf
                        v_micro_lat = -ω_spin_o_dim * y
                        v_micro_roll = ω_spin_o_dim * x + u_heathcote

                        v_lat = v_micro_lat                            # lateral = spin only
                        v_roll = u_slide_o_dim + v_micro_roll          # rolling = macro slide + Heathcote
                        v_mag = sqrt(v_lat^2 + v_roll^2 + 1e-16)

                        mu_local_iso = traction_coefficient(v_mag, trac_A, trac_B, trac_C, trac_D) * ehl_scale_o
                        τ_iso = mu_local_iso * p_local

                        τ_lim = Λ_LSS * p_local
                        τ_lim_eff = sqrt(τ_lim^2 + 1e-24)
                        τ_actual = τ_iso / hypot(1.0, τ_iso / τ_lim_eff) + 1e-6 * τ_iso

                        q_local = τ_actual * v_mag
                        thermal_inertia = 12000.0
                        dT_flash = (1.11 * q_local * sqrt(0.5 * a_o_dim)) /
                                   (thermal_inertia * sqrt(u_roll_eff_o))
                        Na = β_temp * dT_flash
                        thermal_reduction = max(0.25, 1.0 / sqrt(1.0 + Na))

                        τ_final = τ_actual * thermal_reduction

                        dF = τ_final * dx_grid * dy_grid
                        dF_lat = v_mag > 1e-20 ? dF * (v_lat / v_mag) : 0.0
                        dF_roll = v_mag > 1e-20 ? dF * (v_roll / v_mag) : 0.0

                        H_slide_o += dF_roll * u_slide_o_dim    # rolling-direction power
                        H_spin_o += dF_lat * v_micro_lat + dF_roll * v_micro_roll
                    end
                end
            end
        end

        # ── Drag heat [W] ──
        H_drag_visc = F_d_dim * V_ball_dim
        H_drag = H_drag_visc

        # ── EHL Rolling Resistance (Goksem-Hargreaves, matching Python) ──
        churn_rho = ρ_eff
        ν_cst = churn_rho > 0 ? (lub_mu0 / churn_rho) * 1e6 : 10.0
        ndm_factor = 60000.0 / π

        if Q_i_dim > 0 && abs(u_entrain_i) > 1e-6
            # FIX: Goksem-Hargreaves rolling resistance uses the dimensional ball radius
            # as the characteristic length scale, matching Python's dynamics_numba.py.
            R_eff_roll_i = 0.5 * D_ball_dim
            E_ehl = E_prime_dim
            U_i = lub_mu0 * abs(u_entrain_i) / (E_ehl * R_eff_roll_i + 1e-30)
            G_i = lub_alpha_pv * E_ehl
            W_ehl_i = Q_i_dim / (E_ehl * R_eff_roll_i^2 + 1e-30)
            UG_i = U_i * G_i
            if UG_i > 0 && W_ehl_i > 0
                F_roll_i = 4.318 * E_ehl * R_eff_roll_i^2 / G_i * UG_i^0.658 * W_ehl_i^0.0126
                ϕ_ish_i = 1.0 / (1.0 + 1.84e-9 * (abs(u_entrain_i) * ndm_factor)^1.28 * ν_cst^0.64)
                H_drag += F_roll_i * ϕ_ish_i * abs(u_entrain_i)
            end
        end
        if Q_o_dim > 0 && abs(u_entrain_o) > 1e-6
            R_eff_roll_o = 0.5 * D_ball_dim
            E_ehl = E_prime_dim
            U_o = lub_mu0 * abs(u_entrain_o) / (E_ehl * R_eff_roll_o + 1e-30)
            G_o = lub_alpha_pv * E_ehl
            W_ehl_o = Q_o_dim / (E_ehl * R_eff_roll_o^2 + 1e-30)
            UG_o = U_o * G_o
            if UG_o > 0 && W_ehl_o > 0
                F_roll_o = 4.318 * E_ehl * R_eff_roll_o^2 / G_o * UG_o^0.658 * W_ehl_o^0.0126
                ϕ_ish_o = 1.0 / (1.0 + 1.84e-9 * (abs(u_entrain_o) * ndm_factor)^1.28 * ν_cst^0.64)
                H_drag += F_roll_o * ϕ_ish_o * abs(u_entrain_o)
            end
        end

        # ── Churning heat [W] ──
        ω_ball_mag_dim = sqrt(ω[1]^2 + ω[2]^2 + ω[3]^2) * W_scale
        r_ball_dim = D_star / 2 * L_scale
        M_ch_dim = churning_moment(ω_ball_mag_dim, r_ball_dim, ρ_eff, μ_oil)
        H_churn = M_ch_dim * ω_ball_mag_dim

        if j == 1
            # println("  Churning: |ω_body|=$(round(ω_ball_mag_dim, digits=1)), Harris |ω_R|=$(round(abs(ω_R_dim), digits=1)), r_ball=$(round(r_ball_dim*1e3, digits=3))mm, ρ=$(round(ρ_eff, digits=1)), μ=$(round(μ_oil, sigdigits=3))")
            # println("  M_ch=$(round(M_ch_dim*1e6, digits=2))μN·m, H_churn=$(round(H_churn, digits=3))W")
            # H_ehl = H_drag - H_drag_visc
            # println("  H_drag_visc=$(round(H_drag_visc, digits=3))W, H_ehl=$(round(H_ehl, digits=3))W, H_drag_total=$(round(H_drag, digits=3))W")
            # println("  F_d=$(round(F_d_dim, digits=4))N, V_ball=$(round(V_ball_dim, digits=2))m/s")
        end

        # ── Write field output (all in DIMENSIONAL units) ──
        fo[base+FO_Q_I] = Q_i_dim                   # [N]
        fo[base+FO_Q_O] = Q_o_dim                   # [N]
        fo[base+FO_ALPHA_I] = α_i                    # [rad]
        fo[base+FO_ALPHA_O] = α_o                    # [rad]
        fo[base+FO_DELTA_I] = δ_i * L_scale          # [m]
        fo[base+FO_DELTA_O] = δ_o * L_scale          # [m]
        fo[base+FO_A_I] = a_i_dim                    # [m]
        fo[base+FO_B_I] = b_i_dim                    # [m]
        fo[base+FO_A_O] = a_o_dim                    # [m]
        fo[base+FO_B_O] = b_o_dim                    # [m]
        fo[base+FO_U_SLIDE_I] = u_slide_i_dim        # [m/s]
        fo[base+FO_U_SLIDE_O] = u_slide_o_dim        # [m/s]
        fo[base+FO_U_MEAN_I] = u_mean_i * V_scale    # [m/s]
        fo[base+FO_U_MEAN_O] = u_mean_o * V_scale    # [m/s]
        fo[base+FO_H_FILM_I] = h_film_i              # [m]
        fo[base+FO_H_FILM_O] = h_film_o              # [m]
        fo[base+FO_KAPPA_I] = κ_i                    # [-]
        fo[base+FO_KAPPA_O] = κ_o                    # [-]
        fo[base+FO_F_TRAC_I] = F_trac_i_dim          # [N]
        fo[base+FO_F_TRAC_O] = F_trac_o_dim          # [N]
        fo[base+FO_M_SPIN_I] = M_sp_i_dim            # [N·m]
        fo[base+FO_M_SPIN_O] = M_sp_o_dim            # [N·m]
        fo[base+FO_F_DRAG] = F_d_dim                  # [N]
        fo[base+FO_M_CHURN] = M_ch_dim                # [N·m]
        fo[base+FO_OMEGA_X] = ω[1] * W_scale         # [rad/s]
        fo[base+FO_OMEGA_Y] = ω[2] * W_scale         # [rad/s]
        fo[base+FO_OMEGA_Z] = ω[3] * W_scale         # [rad/s]
        fo[base+FO_R_BALL] = r_b * L_scale            # [m]
        fo[base+FO_THETA_DOT] = θ̇_b * W_scale        # [rad/s]
        fo[base+FO_Q0] = real(q)
        v1, v2, v3 = imag_part(q)
        fo[base+FO_Q1] = v1
        fo[base+FO_Q2] = v2
        fo[base+FO_Q3] = v3
        fo[base+FO_H_SLIDE_I] = H_slide_i            # [W]
        fo[base+FO_H_SLIDE_O] = H_slide_o            # [W]
        fo[base+FO_H_SPIN_I] = H_spin_i              # [W]
        fo[base+FO_H_SPIN_O] = H_spin_o              # [W]
        fo[base+FO_H_DRAG] = H_drag                   # [W]
        fo[base+FO_H_CHURN] = H_churn                 # [W]
        fo[base+FO_F_POCKET] = F_pk                   # [N]
        fo[base+FO_W_SPIN_I] = ω_spin_i_dim          # [rad/s]
        fo[base+FO_W_SPIN_O] = ω_spin_o_dim          # [rad/s]
    end

    return fo
end

# ── Helper: Hamrock-Dowson EHL Film Thickness ──
"""
    _film_thickness_hd(u_mean, Q, E_prime, R_eff, a, b, mu_0, alpha_pv, K_th, T_0)

Hamrock-Dowson film thickness with thermal correction.
Returns film thickness h [m].
"""
function _film_thickness_hd(u_mean, Q, E_prime, R_eff, a, b,
    mu_0, alpha_pv, K_th, T_0)
    if Q <= 0 || abs(u_mean) < 1e-15 || R_eff <= 0 || E_prime <= 0
        return 0.0
    end

    # Dimensionless EHL parameters
    U = mu_0 * abs(u_mean) / (E_prime * R_eff)
    G = alpha_pv * E_prime
    W = Q / (E_prime * R_eff^2)

    (U <= 0 || G <= 0 || W <= 0) && return 0.0

    # Ellipticity ratio
    kappa_e = b > 1e-30 ? a / b : 1.0
    ellip_factor = max(1.0 - 0.61 * exp(-0.73 * max(kappa_e, 0.1)), 0.01)

    # Log-domain H-D to avoid underflow
    ln_h_R = log(2.69) + 0.67 * log(U) + 0.53 * log(G) +
             log(ellip_factor) - 0.067 * log(W)

    h_iso = R_eff * exp(ln_h_R)

    # Thermal correction φ_T = 1/(1 + 0.1·Q_m^0.64)
    if K_th > 0 && T_0 > 0
        Q_m = 2.0 * mu_0 * u_mean^2 / (K_th * T_0)
        if Q_m > 1e-10
            phi_T = max(1.0 / (1.0 + 0.1 * Q_m^0.64), 0.01)
        else
            phi_T = 1.0
        end
    else
        phi_T = 1.0
    end

    return max(h_iso * phi_T, 1e-12)
end
