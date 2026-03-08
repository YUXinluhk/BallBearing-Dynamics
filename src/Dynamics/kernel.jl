# =====================================================================
# Dynamics/kernel.jl — ODE RHS for 6-DOF bearing dynamics
#
# Port of dynamics_numba.py::_compute_core (800+ lines).
# Pure function: (du, u, p, t) → nothing.
# Uses @inbounds, @fastmath for performance.
# Quaternion kinematics instead of Euler angles.
# =====================================================================

"""
    ode_rhs!(du, u, params_tuple, t)

In-place ODE right-hand side for the bearing dynamics.

`params_tuple` = (p, ln_re_table, ln_cd_table) where p is the flat parameter vector.
"""
function ode_rhs!(du::AbstractVector, u::AbstractVector,
    params_tuple, t)

    p = params_tuple[1]

    Z = Int(p[P_NBALL])
    D_star = p[P_D]
    d_m_star = p[P_DM]
    f_i = p[P_FI]
    f_o = p[P_FO]
    α₀ = p[P_ALPHA0]
    m_ball = p[P_MBALL]
    J_ball = p[P_JBALL]
    m_ir = p[P_MIR]
    Y_i = p[P_YI]
    Y_o = p[P_YO]
    a_i_star = p[P_AI]
    b_i_star = p[P_BI]
    a_o_star = p[P_AO]
    b_o_star = p[P_BO]
    Σρ_i = p[P_SRI]
    Σρ_o = p[P_SRO]
    D_i_star = p[P_DI]
    D_o_star = p[P_DO]
    μ₀ = p[P_MU0]
    α_pv = p[P_ALPHA_PV]
    trac_A = p[P_TRAC_A]
    trac_B = p[P_TRAC_B]
    trac_C = p[P_TRAC_C]
    trac_D = p[P_TRAC_D]
    μ_spin = p[P_MU_SPIN]
    c_struct = p[P_C_STRUCT]
    ζ = p[P_ZETA]
    ω_ir_star = p[P_OMEGA_IR]
    F_a_star = p[P_FA]
    F_r_star = p[P_FR]
    t_ramp = p[P_T_RAMP]
    ρ_eff = p[P_RHO_EFF]
    μ_oil = p[P_MU_OIL]
    cage_web = p[P_CAGE_WEB]
    R_ball = D_star / 2.0
    drb_i = p[P_DRB_I]
    drb_o = p[P_DRB_O]
    x_gi0_star = p[P_X_GI0]   # inner groove center axial offset ★
    x_go0_star = p[P_X_GO0]   # outer groove center axial offset ★
    # Per-mode damping
    c_ball_t = p[P_C_BALL_TRANS]
    c_ball_o = p[P_C_BALL_ORBIT]
    c_ir_damp = p[P_C_IR_DAMP]
    c_cage_damp = p[P_C_CAGE_DAMP]
    c_tilt_damp = p[P_C_TILT]
    ω_cage_star = p[P_OMEGA_CAGE]
    V_scale = p[P_V_SCALE]    # [m/s] for dimensionalizing velocities
    W_scale = p[P_W_SCALE]    # [rad/s] for dimensionalizing angular velocities
    L_scale = p[P_L_SCALE]    # [m] for dimensionalizing lengths
    Q_scale = p[P_Q_SCALE]    # [N] for dimensionalizing forces

    # Cage params
    cage_mass = p[P_CAGE_MASS]
    cage_Ixx = p[P_CAGE_IXX]
    cage_Iyy = p[P_CAGE_IYY]
    k_pocket = p[P_K_POCKET]
    μ_pocket = p[P_MU_POCKET]
    c_cage = p[P_C_CAGE]

    # ── Load & speed ──
    # Load: full from t=0 (IC is at QS full-load equilibrium)
    F_a = F_a_star
    F_r = F_r_star
    # Speed: C¹ cosine ramp (only speed is ramped)
    ramp = t_ramp > 0 ? 0.5 * (1.0 - cos(π * min(t / t_ramp, 1.0))) : 1.0
    ω_ir = ω_ir_star * ramp

    # ── Read inner race state ──
    ir_p = ir_pos_view(u, Z)
    ir_v = ir_vel_view(u, Z)
    x_ir, y_ir, z_ir = ir_p[1], ir_p[2], ir_p[3]
    γy_ir, γz_ir = ir_p[4], ir_p[5]
    ẋ_ir, ẏ_ir, ż_ir = ir_v[1], ir_v[2], ir_v[3]
    γ̇y_ir, γ̇z_ir = ir_v[4], ir_v[5]

    # ── Read cage state ──
    cg_p = cage_pos_view(u, Z)
    cg_v = cage_vel_view(u, Z)
    x_cg, y_cg, z_cg, θ_cg = cg_p[1], cg_p[2], cg_p[3], cg_p[4]
    ẋ_cg, ẏ_cg, ż_cg, θ̇_cg = cg_v[1], cg_v[2], cg_v[3], cg_v[4]

    # ── Scalar accumulators (zero allocation) ──
    T_u = eltype(u)
    F_ir_x = F_a - c_ir_damp * ẋ_ir
    F_ir_y = zero(T_u) - c_ir_damp * ẏ_ir
    F_ir_z = F_r - c_ir_damp * ż_ir
    M_ir_y = zero(T_u)
    M_ir_z = zero(T_u)

    # Cage: global damping OUTSIDE contact block (C∞ continuous)
    F_cg_x = zero(T_u) - c_cage_damp * ẋ_cg
    F_cg_y = zero(T_u) - c_cage * ẏ_cg
    F_cg_z = zero(T_u) - c_cage * ż_cg
    F_cg_θ = zero(T_u) - c_cage_damp * (θ̇_cg - ω_cage_star)

    N_p = n_pos_dofs(Z)
    N_v = n_vel_dofs(Z)

    # ── Position derivatives = velocities ──
    # Inner race
    @inbounds for i in 1:N_IR_POS
        du[ir_pos_offset()+i-1] = ir_v[i]
    end

    # Cage translational & rotational
    @inbounds for i in 1:N_CAGE_POS
        du[cage_pos_offset(Z)+i-1] = cg_v[i]
    end

    # ── Per-ball loop ──
    @inbounds @fastmath for j in 1:Z
        bp = ball_pos_view(u, j, Z)
        bv = ball_vel_view(u, j, Z)

        x_b, r_b, θ_b = bp[1], bp[2], bp[3]
        q = ball_quat(u, j, Z)
        ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]
        ω_body = ball_omega(u, j, Z)

        sθ, cθ = sincos(θ_b)

        # ── Ball position in inertial frame ──
        y_b = -r_b * sθ
        z_b = r_b * cθ

        # ── Groove center positions (Harris §3.5) ──
        # Strong coordinate coupling: tilt angles γy,γz project into groove center
        R_pitch_i = D_i_star / 2.0
        r_gc_i = R_pitch_i + y_ir * (-sθ) + z_ir * cθ  # + IR lateral displacement
        # Axial: x_gi0★ + IR axial + tilt projection (breakthrough 1: bend-twist coupling)
        x_gc_i = x_gi0_star + x_ir + R_pitch_i * (cθ * γy_ir + sθ * γz_ir)

        r_gc_o = D_o_star / 2   # outer race is fixed
        x_gc_o = x_go0_star     # outer groove axial offset (non-zero!)

        # ── Contact geometry — 3D normal vector approach (matches Python) ──
        # Vectors from groove center to ball center
        dx_i = x_gc_i - x_b     # axial gap inner (x_ir only in x_gc_i, not double-counted)
        dr_i = r_gc_i - r_b     # radial gap inner

        dx_o = x_gc_o - x_b     # axial gap outer (includes x_go offset)
        dr_o = r_gc_o - r_b     # radial gap outer

        # Distance from ball center to groove center
        L_i = sqrt(dx_i^2 + dr_i^2 + 1e-30)
        L_o = sqrt(dx_o^2 + dr_o^2 + 1e-30)

        # Unit normal vectors: groove center → ball center direction
        inv_Li = 1.0 / L_i
        inv_Lo = 1.0 / L_o
        nfi_x = dx_i * inv_Li    # = sin(α_i)
        nfi_r = dr_i * inv_Li    # = cos(α_i)
        nfo_x = dx_o * inv_Lo    # = sin(α_o) (negative for outer)
        nfo_r = dr_o * inv_Lo    # = cos(α_o) (negative for outer)

        # Penetration (elastic deformation) — SAME formula for both
        δ_i = L_i - drb_i    # positive = contact
        δ_o = L_o - drb_o    # positive = contact (FIX: was drb_o - L_o)

        # Smooth Hertz (C∞ continuous)
        δ_i_sm = smooth_hertz_delta(δ_i)
        δ_o_sm = smooth_hertz_delta(δ_o)

        # ── Hertz contact loads (x*sqrt(x) ≡ x^1.5 but 3× faster) ──
        Q_i = Y_i * δ_i_sm * sqrt(δ_i_sm)
        Q_o = Y_o * δ_o_sm * sqrt(δ_o_sm)

        # Contact angles (for traction/spin, derived from normals)
        sαi = nfi_x    # sin(α_i)
        cαi = nfi_r    # cos(α_i)
        sαo = nfo_x    # sin(α_o)
        cαo = nfo_r    # cos(α_o)

        # ── Contact ellipse dimensions (cbrt ≡ ^(1/3) but native HW) ──
        c_i = Σρ_i > 1e-20 ? cbrt(3Q_i / (2 * p[P_E_PRIME] * Σρ_i)) : zero(T_u)
        a_i = a_i_star * c_i
        b_i = b_i_star * c_i
        c_o = Σρ_o > 1e-20 ? cbrt(3Q_o / (2 * p[P_E_PRIME] * Σρ_o)) : zero(T_u)
        a_o = a_o_star * c_o
        b_o = b_o_star * c_o

        # ── Kinematics: surface velocities (nondim) ──
        R_ball_k = D_star / 2  # ball radius (nondim)
        # Inner race surface velocity at contact point (not groove center!)
        v_ir_θ = (r_b - R_ball_k * cαi) * ω_ir  # contact point radius (nondim)
        v_ball_θ = r_b * θ̇_b            # ball orbital velocity (nondim)

        # Ball surface velocity at contact point (includes ball rotation!)
        # ω_roll = ωx·cos(α) - ωz·sin(α) — rolling component of ball angular velocity
        ω_roll_i = ω_body[1] * cαi - ω_body[3] * sαi
        ω_roll_o = ω_body[1] * cαo - ω_body[3] * sαo

        # Ball surface velocity in tangential direction at contact:
        # Left-hand (x,θ,r) system: ê_x × ê_θ = -ê_r ⟹ v_surf_θ = -R·ω_roll
        v_ball_surface_i = v_ball_θ - R_ball_k * ω_roll_i  # nondim
        v_ball_surface_o = v_ball_θ - R_ball_k * ω_roll_o  # nondim

        # Dimensionalize for traction coefficient (expects m/s)
        v_ir_θ_dim = v_ir_θ * V_scale
        v_ball_θ_dim = v_ball_θ * V_scale

        # Slide velocity: creep = race surface - ball contact surface
        u_slide_i_dim = abs(v_ir_θ_dim - v_ball_surface_i * V_scale)
        u_slide_o_dim = abs(v_ball_surface_o * V_scale)  # outer race is fixed

        # Spin velocity (dimensional rad/s)
        ω_ir_dim = ω_ir * W_scale
        ω_race_n_i = ω_ir_dim * sαi
        ω_ball_n_i = (ω_body[1] * sαi + ω_body[3] * cαi) * W_scale  # ω⃗·n̂ = ωₓsinα + ω_r·cosα
        ω_spin_i = ω_race_n_i - ω_ball_n_i

        ω_race_n_o = zero(T_u)  # outer race is fixed
        ω_ball_n_o = (ω_body[1] * sαo + ω_body[3] * cαo) * W_scale  # dot product (always additive)
        ω_spin_o = ω_race_n_o - ω_ball_n_o

        # ── TEHD traction forces (Harris paradox killer) ──
        # Dimensionalize for TEHD model
        Q_i_dim = Q_i * Q_scale
        Q_o_dim = Q_o * Q_scale
        a_i_dim_t = a_i * L_scale
        b_i_dim_t = b_i * L_scale
        a_o_dim_t = a_o * L_scale
        b_o_dim_t = b_o * L_scale

        # Mean contact pressure P = Q / (π·a·b)
        P_mean_i = Q_i_dim / (π * a_i_dim_t * b_i_dim_t + 1e-16)
        P_mean_o = Q_o_dim / (π * a_o_dim_t * b_o_dim_t + 1e-16)

        # Entrainment velocities
        v_ball_θ_dim = v_ball_θ * V_scale
        u_entrain_i = 0.5 * (v_ir_θ_dim + v_ball_θ_dim)
        u_entrain_o = 0.5 * v_ball_θ_dim

        # TEHD traction coefficient (LSS + flash temp + wave speed limit)
        Lambda_LSS = p[P_LAMBDA_LSS]
        beta_temp = p[P_BETA_TEMP]
        κ_i = tehd_traction_coefficient(u_slide_i_dim, P_mean_i, a_i_dim_t, u_entrain_i,
            trac_A, trac_B, trac_C, trac_D, Lambda_LSS, beta_temp)
        κ_o = tehd_traction_coefficient(u_slide_o_dim, P_mean_o, a_o_dim_t, u_entrain_o,
            trac_A, trac_B, trac_C, trac_D, Lambda_LSS, beta_temp)

        F_trac_i = (κ_i * Q_i_dim) / Q_scale   # nondim
        F_trac_o = (κ_o * Q_o_dim) / Q_scale   # nondim

        # Traction direction — C∞ smooth (tanh replaces sign)
        v_slip_i = v_ir_θ_dim - v_ball_surface_i * V_scale
        v_slip_o = -v_ball_surface_o * V_scale
        sign_slide_i = tanh(v_slip_i / 0.01)
        sign_slide_o = tanh(v_slip_o / 0.01)

        # ── Spin moment ──
        # Contact ellipse in dimensional [m] for spin moment
        a_i_dim = a_i * L_scale
        b_i_dim = b_i * L_scale
        a_o_dim = a_o * L_scale
        b_o_dim = b_o * L_scale
        M_spin_i = spin_moment(ω_spin_i, Q_i_dim, a_i_dim, b_i_dim, μ_spin, p[P_E2_I]) / (Q_scale * L_scale)
        M_spin_o = spin_moment(ω_spin_o, Q_o_dim, a_o_dim, b_o_dim, μ_spin, p[P_E2_O]) / (Q_scale * L_scale)

        # ── Drag force ──
        # Drag: dimensionalize velocity and ball diameter, then nondim the result
        V_ball_dim = sqrt(ẋ_b^2 + (r_b * θ̇_b)^2 + ṙ_b^2) * V_scale
        D_ball_dim = D_star * L_scale
        F_drag = drag_force(V_ball_dim, D_ball_dim, ρ_eff, μ_oil, cage_web * L_scale) / Q_scale

        # ── Ball force balance ──
        # Normal contact forces: Q × unit_normal (FIX: no sign games, direction from normal)
        F_n_ix = Q_i * nfi_x    # inner normal, axial component
        F_n_ir = Q_i * nfi_r    # inner normal, radial component
        F_n_ox = Q_o * nfo_x    # outer normal, axial (FIX: was -Q_o * sαo)
        F_n_or = Q_o * nfo_r    # outer normal, radial (FIX: was -Q_o * cαo)

        # Traction forces (tangential, in θ direction)
        F_trac_tang_i = F_trac_i * sign_slide_i
        F_trac_tang_o = F_trac_o * sign_slide_o

        # Total ball forces in cylindrical coords (damping applied in acc equations)
        F_ball_x = F_n_ix + F_n_ox
        F_ball_r = F_n_ir + F_n_or
        F_ball_θ = F_trac_tang_i + F_trac_tang_o - F_drag * tanh(r_b * θ̇_b * V_scale / 0.01)

        # Ball/cage pocket interaction (simplified)
        pocket_θ = θ_cg + (j - 1) * 2π / Z
        Δθ_pocket = θ_b - pocket_θ
        Δθ_pocket -= 2π * round(Δθ_pocket / (2π))
        gap_pocket = r_b * Δθ_pocket
        pocket_clr = p[P_POCKET_CLR]
        pen_pocket = abs(gap_pocket) - pocket_clr
        pen_pocket_sm = smooth_hertz_delta(pen_pocket)
        F_pocket = k_pocket * pen_pocket_sm * sqrt(pen_pocket_sm)
        sign_pocket = tanh(gap_pocket / 1e-6)
        F_ball_θ -= sign_pocket * F_pocket
        F_cg_θ += sign_pocket * F_pocket * r_b  # cage torque

        # (c_lin damping removed — already applied via c_ball_t at acceleration level)

        # ── 3D Rolling torque: r_contact × F_traction ──
        # In Julia's 2D frame (0=axial, 1=tangent, 2=radial):
        #   normal = (nfi_x, 0, nfi_r), traction = (0, F_tang, 0)
        #   contact_pt = -R_ball * normal = (-R·sα, 0, -R·cα)
        #   τ = r × F = (R·cα·F, 0, -R·sα·F) (dimensional → nondim)
        R_ball_dim = D_star * L_scale / 2
        moment_scale_inv = 1.0 / (Q_scale * L_scale)  # nondim moment scale

        # Inner contact traction torque (dimensional cross product → nondim)
        F_tang_i_dim = F_trac_tang_i * Q_scale  # dimensional [N]
        τ_roll_i_x = -R_ball_dim * cαi * F_tang_i_dim * moment_scale_inv   # = -R·cos(α_i)·F / (Q·L)
        τ_roll_i_z = R_ball_dim * sαi * F_tang_i_dim * moment_scale_inv    # = +R·sin(α_i)·F / (Q·L)

        # Outer contact traction torque
        F_tang_o_dim = F_trac_tang_o * Q_scale  # dimensional [N]
        τ_roll_o_x = -R_ball_dim * cαo * F_tang_o_dim * moment_scale_inv
        τ_roll_o_z = R_ball_dim * sαo * F_tang_o_dim * moment_scale_inv

        # Ball moment: rolling torque + spin moment
        M_ball_x = τ_roll_i_x + τ_roll_o_x + M_spin_i * sαi + M_spin_o * sαo
        M_ball_y = zero(T_u)
        M_ball_z = τ_roll_i_z + τ_roll_o_z + M_spin_i * cαi + M_spin_o * cαo  # spin moment r-component = +M·cosα

        # Rayleigh spin damping (from params.jl per-mode model)
        c_ball_spin = p[P_C_BALL_SPIN]
        M_ball_x -= c_ball_spin * ω_body[1]
        M_ball_y -= c_ball_spin * ω_body[2]
        M_ball_z -= c_ball_spin * ω_body[3]

        # Euler transport: -(Ω_orbit × H) = (0, I·θ̇·ω_z, -I·θ̇·ω_y)
        M_ball_y -= J_ball * θ̇_b * ω_body[3]    # y ↔ ê_θ (tangential)
        M_ball_z += J_ball * θ̇_b * ω_body[2]    # z ↔ ê_r (radial)
        # ── Write ball derivatives ──
        # Position derivatives
        off_p = ball_pos_offset(j)
        du[off_p] = ẋ_b
        du[off_p+1] = ṙ_b
        du[off_p+2] = θ̇_b

        # Quaternion derivative with Baumgarte stabilization (replaces DiscreteCallback)
        # q = (s, v1, v2, v3) stored at u[off_p+3 .. off_p+6]
        q_s = u[off_p+3]
        q_1 = u[off_p+4]
        q_2 = u[off_p+5]
        q_3 = u[off_p+6]
        q_norm_sq = q_s^2 + q_1^2 + q_2^2 + q_3^2
        λ_q = 100.0 * (1.0 - q_norm_sq)   # Baumgarte attractor → |q|=1
        q_dot = quat_derivative(q, ω_body)
        du[off_p+3] = real(q_dot) + λ_q * q_s
        v1, v2, v3 = imag_part(q_dot)
        du[off_p+4] = v1 + λ_q * q_1
        du[off_p+5] = v2 + λ_q * q_2
        du[off_p+6] = v3 + λ_q * q_3

        # Velocity derivatives (accelerations) — per-mode damping applied here
        off_v = ball_vel_offset(j, Z)
        du[off_v] = F_ball_x / m_ball - c_ball_t * ẋ_b
        du[off_v+1] = F_ball_r / m_ball + r_b * θ̇_b^2 - c_ball_t * ṙ_b
        du[off_v+2] = (F_ball_θ / m_ball - 2 * ṙ_b * θ̇_b) / (r_b + 1e-30) - c_ball_o * (θ̇_b - ω_cage_star)
        du[off_v+3] = M_ball_x / J_ball
        du[off_v+4] = M_ball_y / J_ball
        du[off_v+5] = M_ball_z / J_ball

        # ── Accumulate on inner race (Newton's 3rd law) ──
        F_ir_x -= F_n_ix
        F_ir_y -= (-F_n_ir * sθ - F_trac_tang_i * cθ)
        F_ir_z -= (F_n_ir * cθ - F_trac_tang_i * sθ)

        # ── IR tilt moments from contact forces ──
        M_ir_y += -F_n_ix * (D_i_star / 2) * cθ + x_gc_i * (F_n_ir * cθ - F_trac_tang_i * sθ)
        M_ir_z += -F_n_ix * (D_i_star / 2) * sθ + x_gc_i * (F_n_ir * sθ + F_trac_tang_i * cθ)
    end

    # ── Inner race acceleration ──
    off_irv = ir_vel_offset(Z)
    du[off_irv] = F_ir_x / m_ir
    du[off_irv+1] = F_ir_y / m_ir
    du[off_irv+2] = F_ir_z / m_ir
    # IR tilt: moment from contact forces + tilt damping
    I_ir_tilt = m_ir * (D_i_star / 2)^2   # approximate tilt inertia
    du[off_irv+3] = (M_ir_y - c_tilt_damp * γ̇y_ir) / (I_ir_tilt + 1e-30)
    du[off_irv+4] = (M_ir_z - c_tilt_damp * γ̇z_ir) / (I_ir_tilt + 1e-30)

    # ── Cage pilot force (C∞ smooth) ──
    ecc = sqrt(y_cg^2 + z_cg^2)
    pilot_clr = p[P_PILOT_CLR]
    δ_pilot = smooth_hertz_delta(ecc - pilot_clr)
    F_pilot = p[P_K_PILOT] * δ_pilot * sqrt(δ_pilot)
    ecc_safe = ecc + 1e-30
    ey, ez = y_cg / ecc_safe, z_cg / ecc_safe
    F_cg_y -= F_pilot * ey
    F_cg_z -= F_pilot * ez

    # Pilot drag torque (smooth sign)
    Δω_pilot = θ̇_cg - (p[P_PILOT_IR] > 0.5 ? ω_ir : zero(T_u))
    R_pilot = p[P_PILOT_IR] > 0.5 ? p[P_CAGE_IR] : p[P_CAGE_OR]
    F_cg_θ -= p[P_MU_PILOT] * F_pilot * tanh(Δω_pilot / 1e-3) * R_pilot

    # Cage acceleration
    off_cv = cage_vel_offset(Z)
    du[off_cv] = F_cg_x / cage_mass
    du[off_cv+1] = F_cg_y / cage_mass
    du[off_cv+2] = F_cg_z / cage_mass
    du[off_cv+3] = F_cg_θ / cage_Ixx

    return nothing
end
