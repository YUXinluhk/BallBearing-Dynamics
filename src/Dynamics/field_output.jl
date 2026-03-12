# =====================================================================
# Dynamics/field_output.jl — Field output extraction (two-pass approach)
# 完全镜像同步 kernel 动力学内核 (deepthink.md 参考架构)
# =====================================================================

# Field output indices per ball
const N_FIELD_PER_BALL = 42
const FO_Q_I = 1;
const FO_Q_O = 2;
const FO_ALPHA_I = 3;
const FO_ALPHA_O = 4;
const FO_DELTA_I = 5;
const FO_DELTA_O = 6;
const FO_A_I = 7;
const FO_B_I = 8;
const FO_A_O = 9;
const FO_B_O = 10;
const FO_U_SLIDE_I = 11;
const FO_U_SLIDE_O = 12;
const FO_U_MEAN_I = 13;
const FO_U_MEAN_O = 14;
const FO_H_FILM_I = 15;
const FO_H_FILM_O = 16;
const FO_KAPPA_I = 17;
const FO_KAPPA_O = 18;
const FO_F_TRAC_I = 19;
const FO_F_TRAC_O = 20;
const FO_M_SPIN_I = 21;
const FO_M_SPIN_O = 22;
const FO_F_DRAG = 23;
const FO_M_CHURN = 24;
const FO_OMEGA_X = 25;
const FO_OMEGA_Y = 26;
const FO_OMEGA_Z = 27;
const FO_R_BALL = 28;
const FO_THETA_DOT = 29;
const FO_Q0 = 30;
const FO_Q1 = 31;
const FO_Q2 = 32;
const FO_Q3 = 33;
const FO_H_SLIDE_I = 34;
const FO_H_SLIDE_O = 35;
const FO_H_SPIN_I = 36;
const FO_H_SPIN_O = 37;
const FO_H_DRAG = 38;
const FO_H_CHURN = 39;
const FO_F_POCKET = 40;
const FO_W_SPIN_I = 41;
const FO_W_SPIN_O = 42;

"""
    field_output_kernel(t, u, p) → Vector{Float64}

Extract diagnostic field outputs from state vector u at time t.
完全镜像同步：与 kernel.jl 使用完全相同的接触几何和运动学公式。
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
    ρ_eff = p[P_RHO_EFF]
    μ_oil = p[P_MU_OIL]
    cage_web = p[P_CAGE_WEB]
    D_i_star = p[P_DI]
    D_o_star = p[P_DO]
    x_gi0_star = p[P_X_GI0]
    x_go0_star = p[P_X_GO0]
    V_scale = p[P_V_SCALE]
    W_scale = p[P_W_SCALE]
    L_scale = p[P_L_SCALE]
    Q_scale = p[P_Q_SCALE]

    ramp = p[P_T_RAMP] > 0 ? 0.5 * (1.0 - cos(π * min(t / p[P_T_RAMP], 1.0))) : 1.0
    ω_ir = p[P_OMEGA_IR] * ramp
    composite_E_prime = p[P_E_PRIME] * Q_scale / L_scale^2

    lub_mu0 = p[P_MU0]
    lub_alpha_pv = p[P_ALPHA_PV]
    lub_K_th = p[P_KTH]
    lub_T0 = p[P_T0]
    Λ_LSS = p[P_LAMBDA_LSS]
    β_temp = p[P_BETA_TEMP]
    f_i_geom = p[P_FI]
    f_o_geom = p[P_FO]

    # ── 内外圈状态 (完全镜像 kernel.jl) ──
    ir_p = ir_pos_view(u, Z)
    ir_v = ir_vel_view(u, Z)
    x_ir, y_ir, z_ir, γy_ir, γz_ir = ir_p[1], ir_p[2], ir_p[3], ir_p[4], ir_p[5]
    ẋ_ir, ẏ_ir, ż_ir, γ̇y_ir, γ̇z_ir = ir_v[1], ir_v[2], ir_v[3], ir_v[4], ir_v[5]

    cg_v = cage_vel_view(u, Z)
    θ̇_cg = cg_v[4]
    cg_p = cage_pos_view(u, Z)
    θ_cg = cg_p[4]

    # Thermal node temperatures (mirroring kernel.jl)
    tv = thermal_view(u, Z)
    T_ir_t, T_or_t = tv[1], tv[2]
    T_ball_t = tv[3]
    T_oil_t = tv[4]

    # Phase 4b: Ball thermal expansion (mirroring kernel.jl)
    R_ball_k_nom = D_star / 2.0
    δr_ball = p[P_CTE] * (T_ball_t - p[P_T_REF]) * R_ball_k_nom
    R_ball_k_global = R_ball_k_nom + δr_ball
    # Dynamic viscosity (bidirectional Roelands, mirroring kernel.jl)
    μ_oil_dyn = μ_oil * exp(-β_temp * (T_oil_t - lub_T0))

    @inbounds for j in 1:Z
        base = (j - 1) * N_FIELD_PER_BALL
        bp = ball_pos_view(u, j, Z)
        bv = ball_vel_view(u, j, Z)
        x_b, r_b, θ_b = bp[1], bp[2], bp[3]
        ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]
        ω_orb = ball_omega(u, j, Z)
        q = ball_quat(u, j, Z)
        sθ, cθ = sincos(θ_b)

        ω_orb_x, ω_orb_r, ω_orb_θ = ω_orb[1], ω_orb[2], ω_orb[3]

        # 【核心修复1】完整接触几何：内圈倾覆投影同步 (与 kernel.jl 完全镜像)
        # Phase 4: 热膨胀动态游隙修正 — 参考温度为初始运行温度 T0 (镜像 kernel.jl)
        δr_ir = p[P_CTE] * (T_ir_t - p[P_T_REF]) * (D_i_star / 2.0)
        δr_or = p[P_CTE] * (T_or_t - p[P_T_REF]) * (D_o_star / 2.0)
        R_pitch_i = D_i_star / 2.0 + δr_ir
        tilt_proj = cθ * γy_ir + sθ * γz_ir
        x_gc_i = x_gi0_star + x_ir + R_pitch_i * tilt_proj
        r_gc_i = R_pitch_i + y_ir * (-sθ) + z_ir * cθ - x_gi0_star * tilt_proj
        x_gc_o = x_go0_star
        r_gc_o = D_o_star / 2.0 + δr_or

        dx_i = x_gc_i - x_b
        dr_i = r_gc_i - r_b
        dx_o = x_gc_o - x_b
        dr_o = r_gc_o - r_b
        L_i = sqrt(dx_i^2 + dr_i^2 + 1e-30)
        L_o = sqrt(dx_o^2 + dr_o^2 + 1e-30)

        nfi_x, nfi_r = dx_i / L_i, dr_i / L_i
        nfo_x, nfo_r = dx_o / L_o, dr_o / L_o

        # Smooth Hertz penetration (同 kernel.jl smooth_hertz_delta 等价连续形式)
        δ_i_sm = 0.5 * (L_i - drb_i + sqrt((L_i - drb_i)^2 + 1e-24))
        δ_o_sm = 0.5 * (L_o - drb_o + sqrt((L_o - drb_o)^2 + 1e-24))
        Q_i = Y_i * δ_i_sm * sqrt(δ_i_sm)
        Q_o = Y_o * δ_o_sm * sqrt(δ_o_sm)

        α_i = atan(dx_i, dr_i + 1e-30)
        α_o = atan(dx_o, -(dr_o) + 1e-30)
        c_i = Σρ_i > 1e-20 ? cbrt(3Q_i / (2 * p[P_E_PRIME] * Σρ_i)) : 0.0
        a_i_dim, b_i_dim = a_i_s * c_i * L_scale, b_i_s * c_i * L_scale
        c_o = Σρ_o > 1e-20 ? cbrt(3Q_o / (2 * p[P_E_PRIME] * Σρ_o)) : 0.0
        a_o_dim, b_o_dim = a_o_s * c_o * L_scale, b_o_s * c_o * L_scale

        R_ball_k = R_ball_k_global  # includes ball thermal expansion
        r_contact_i = r_b - R_ball_k * nfi_r
        x_contact_i = x_b - R_ball_k * nfi_x

        # 【核心修复2】与 kernel.jl 完全一致的接触点线速度公式
        # 内圈在接触点的动态挤压线速度 (包含倾覆振动产生的速度)
        v_ir_θ = r_contact_i * ω_ir + γ̇y_ir * (x_contact_i - x_ir) * sθ - γ̇z_ir * (x_contact_i - x_ir) * cθ - ẏ_ir * cθ - ż_ir * sθ
        v_ball_θ = r_b * θ̇_b
        u_mean_i = 0.5 * (v_ir_θ + v_ball_θ) + 1e-30
        u_mean_o = 0.5 * v_ball_θ + 1e-30

        R_ball_dim = R_ball_k * L_scale
        ω_ir_dim_val = ω_ir * W_scale

        # 【核心修复3】使用 ODE 的 ω_orb (非 Harris 混合模型)，与 kernel.jl 完全一致
        ω_roll_i = ω_orb_x * nfi_r - ω_orb_r * nfi_x
        ω_roll_o = ω_orb_x * nfo_r - ω_orb_r * nfo_x

        v_ball_surface_i = v_ball_θ - R_ball_k * ω_roll_i
        v_ball_surface_o = v_ball_θ - R_ball_k * ω_roll_o

        v_ir_θ_dim = v_ir_θ * V_scale
        # 【Bug8修复】保留有符号滑移速度，维持伽利略对称性
        u_slide_i_dim = v_ir_θ_dim - v_ball_surface_i * V_scale
        u_slide_o_dim = -v_ball_surface_o * V_scale

        # 自旋速度投影 (同 kernel.jl)
        ω_spin_i_dim = ω_ir_dim_val * nfi_x - (ω_orb_x * nfi_x + ω_orb_r * nfi_r) * W_scale
        ω_spin_o_dim = 0.0 - (ω_orb_x * nfo_x + ω_orb_r * nfo_r) * W_scale

        v_surface_i_θ = v_ball_surface_i * V_scale
        v_surface_o_θ = v_ball_surface_o * V_scale
        Q_i_dim = Q_i * Q_scale
        Q_o_dim = Q_o * Q_scale

        E_prime_dim = 2.0 * composite_E_prime  # Hamrock-Dowson: E' = 2E*

        # 【Bug1修复】EHL卷吸速度：接触点共转参考系
        r_contact_o_fo = r_b - R_ball_k * nfo_r
        v_cp_i_dim_fo = r_contact_i * θ̇_b * V_scale
        v_cp_o_dim_fo = r_contact_o_fo * θ̇_b * V_scale
        u_entrain_i = abs(0.5 * (v_ir_θ_dim + v_surface_i_θ) - v_cp_i_dim_fo)
        u_entrain_o = abs(0.5 * (0.0 + v_surface_o_θ) - v_cp_o_dim_fo)

        # ── Hamrock-Dowson Film Thickness ──
        R_eff = R_ball_dim
        # 【Bug3修复】使用动态油温和热软化粘度
        h_film_i = _film_thickness_hd(u_entrain_i, Q_i_dim, E_prime_dim, R_eff,
            (b_i_dim > 0 ? a_i_dim / b_i_dim : 1.0), μ_oil_dyn, lub_alpha_pv, lub_K_th, T_oil_t)
        h_film_o = _film_thickness_hd(u_entrain_o, Q_o_dim, E_prime_dim, R_eff,
            (b_o_dim > 0 ? a_o_dim / b_o_dim : 1.0), μ_oil_dyn, lub_alpha_pv, lub_K_th, T_oil_t)

        # 【核心修复4】使用 TEHD 牵引系数 (完全镜像 kernel.jl，非旧的基础 traction_coefficient)
        P_mean_i = Q_i_dim > 0 ? Q_i_dim / (π * a_i_dim * b_i_dim + 1e-16) : 0.0
        P_mean_o = Q_o_dim > 0 ? Q_o_dim / (π * a_o_dim * b_o_dim + 1e-16) : 0.0
        κ_i = tehd_traction_coefficient(u_slide_i_dim, P_mean_i, b_i_dim, u_entrain_i, trac_A, trac_B, trac_C, trac_D, Λ_LSS, β_temp)
        κ_o = tehd_traction_coefficient(u_slide_o_dim, P_mean_o, b_o_dim, u_entrain_o, trac_A, trac_B, trac_C, trac_D, Λ_LSS, β_temp)

        v_roll_nominal = abs(ω_ir_dim_val) * d_m_star * L_scale * 0.5
        entrain_factor = v_roll_nominal^2 / (v_roll_nominal^2 + 0.01)
        κ_i *= entrain_factor
        κ_o *= entrain_factor

        # (Stribeck 分载已移除，与 kernel.jl 完全镜像)

        F_trac_i_dim = κ_i * Q_i_dim
        F_trac_o_dim = κ_o * Q_o_dim

        M_sp_i_dim = spin_moment(ω_spin_i_dim, Q_i_dim, a_i_dim, b_i_dim, μ_spin, p[P_E2_I])
        M_sp_o_dim = spin_moment(ω_spin_o_dim, Q_o_dim, a_o_dim, b_o_dim, μ_spin, p[P_E2_O])

        # ── 3D 流体拖曳力 ──
        V_ball_dim = sqrt(ẋ_b^2 + (r_b * θ̇_b)^2 + ṙ_b^2) * V_scale
        D_ball_dim = D_star * L_scale
        F_d_dim = drag_force(V_ball_dim, D_ball_dim, ρ_eff, μ_oil_dyn, cage_web * L_scale)

        # ── 兜孔接触力 ──
        pocket_θ = θ_cg + (j - 1) * 2π / Z
        Δθ = θ_b - pocket_θ
        Δθ -= 2π * round(Δθ / (2π))
        gap = r_b * Δθ
        pen = abs(gap) - p[P_POCKET_CLR]
        F_pk = pen > 0 ? p[P_K_POCKET] * pen * Q_scale : 0.0

        # =============================================================
        # 2D 接触椭圆热积分 (11×11 网格)
        # Archard 闪温 + Nahme TEHD 热软化 + Bair-Winer 极限剪应力
        # =============================================================
        H_slide_i, H_slide_o = 0.0, 0.0
        H_spin_i, H_spin_o = 0.0, 0.0

        # [V6] 7-point Gauss-Legendre (u in [0,1]) and 12-point Fourier (θ in [0,2π]) mappings
        nodes_u = (0.02544604382862075, 0.12923440720030275, 0.29707742431130140,
            0.5, 0.70292257568869860, 0.87076559279969725, 0.97455395617137925)
        weights_u = (0.06474248308443485, 0.13985269574463830, 0.19091502525255945,
            0.20897959183673470, 0.19091502525255945, 0.13985269574463830, 0.06474248308443485)
        N_theta = 12
        dtheta = 2.0 * π / N_theta

        ω_roll_abs = sqrt(ω_orb_x^2 + ω_orb_r^2 + ω_orb_θ^2 + 1e-16) * W_scale

        # --- 内圈接触 84-点 积分 (V6) ---
        if a_i_dim > 0 && b_i_dim > 0 && Q_i_dim > 0
            p_max_i = 3.0 * Q_i_dim / (2.0 * π * a_i_dim * b_i_dim)
            # 【Bug2修复】Heathcote共形系数：乘以 (2f-1)/(2f) 而非除以 (2f-1)
            conf_factor_i = (2.0 * f_i_geom - 1.0) / (2.0 * f_i_geom)
            u_roll_eff_i = sqrt(u_entrain_i^2 + 1e-6)
            film_atten_i = 0.7 + 0.3 * (Q_i_dim > 0 && h_film_i > 0.0 ? exp(-0.6 * (h_film_i / 1e-7)^1.5) : 1.0)

            for i_u in 1:7
                u_nd = nodes_u[i_u]
                wu = weights_u[i_u]
                ρ_val = sqrt(1.0 - u_nd^2)
                p_local = p_max_i * u_nd
                τ_lim_eff = sqrt((Λ_LSS * p_local)^2 + 1e-24)

                for i_theta in 1:N_theta
                    θ = (i_theta - 0.5) * dtheta
                    s_th, c_th = sincos(θ)
                    x_c = a_i_dim * ρ_val * c_th
                    y_c = b_i_dim * ρ_val * s_th

                    u_heathcote = ω_roll_abs * (x_c^2 / (2.0 * R_ball_dim)) * conf_factor_i  # 【Bug2修复】
                    v_micro_lat = -ω_spin_i_dim * y_c
                    v_micro_roll = ω_spin_i_dim * x_c + u_heathcote

                    v_lat = v_micro_lat
                    v_roll = u_slide_i_dim + v_micro_roll  # 【Bug8修复】u_slide_i_dim已为有符号
                    v_mag = sqrt(v_lat^2 + v_roll^2 + 1e-16)

                    mu_local_iso = traction_coefficient(v_mag, trac_A, trac_B, trac_C, trac_D) * film_atten_i
                    τ_iso = mu_local_iso * p_local
                    τ_actual = τ_iso / hypot(1.0, τ_iso / τ_lim_eff)

                    q_local = τ_actual * v_mag
                    dT_flash = (1.11 * q_local * sqrt(b_i_dim)) / (12000.0 * sqrt(u_roll_eff_i))
                    thermal_reduction = max(0.25, 1.0 / sqrt(1.0 + β_temp * dT_flash))

                    τ_final = τ_actual * thermal_reduction
                    dF = τ_final * (a_i_dim * b_i_dim * u_nd * wu * dtheta)
                    dF_lat = v_mag > 1e-20 ? dF * (v_lat / v_mag) : 0.0
                    dF_roll = v_mag > 1e-20 ? dF * (v_roll / v_mag) : 0.0

                    H_slide_i += dF_roll * u_slide_i_dim  # 【Bug8修复】有符号自然积分
                    H_spin_i += dF_lat * v_micro_lat + dF_roll * v_micro_roll
                end
            end
        end

        # --- 外圈接触 84-点 积分 (V6) ---
        if a_o_dim > 0 && b_o_dim > 0 && Q_o_dim > 0
            p_max_o = 3.0 * Q_o_dim / (2.0 * π * a_o_dim * b_o_dim)
            # 【Bug2修复】Heathcote共形系数：乘以 (2f-1)/(2f) 而非除以 (2f-1)
            conf_factor_o = (2.0 * f_o_geom - 1.0) / (2.0 * f_o_geom)
            u_roll_eff_o = sqrt(u_entrain_o^2 + 1e-6)
            film_atten_o = 0.7 + 0.3 * (Q_o_dim > 0 && h_film_o > 0.0 ? exp(-0.6 * (h_film_o / 1e-7)^1.5) : 1.0)

            for i_u in 1:7
                u_nd = nodes_u[i_u]
                wu = weights_u[i_u]
                ρ_val = sqrt(1.0 - u_nd^2)
                p_local = p_max_o * u_nd
                τ_lim_eff = sqrt((Λ_LSS * p_local)^2 + 1e-24)

                for i_theta in 1:N_theta
                    θ = (i_theta - 0.5) * dtheta
                    s_th, c_th = sincos(θ)
                    x_c = a_o_dim * ρ_val * c_th
                    y_c = b_o_dim * ρ_val * s_th

                    u_heathcote = ω_roll_abs * (x_c^2 / (2.0 * R_ball_dim)) * conf_factor_o  # 【Bug2修复】
                    v_micro_lat = -ω_spin_o_dim * y_c
                    v_micro_roll = ω_spin_o_dim * x_c + u_heathcote

                    v_lat = v_micro_lat
                    v_roll = u_slide_o_dim + v_micro_roll  # 【Bug8修复】u_slide_o_dim已为有符号
                    v_mag = sqrt(v_lat^2 + v_roll^2 + 1e-16)

                    mu_local_iso = traction_coefficient(v_mag, trac_A, trac_B, trac_C, trac_D) * film_atten_o
                    τ_iso = mu_local_iso * p_local
                    τ_actual = τ_iso / hypot(1.0, τ_iso / τ_lim_eff)

                    q_local = τ_actual * v_mag
                    dT_flash = (1.11 * q_local * sqrt(0.5 * a_o_dim)) / (12000.0 * sqrt(u_roll_eff_o))
                    thermal_reduction = max(0.25, 1.0 / sqrt(1.0 + β_temp * dT_flash))

                    τ_final = τ_actual * thermal_reduction
                    dF = τ_final * (a_o_dim * b_o_dim * u_nd * wu * dtheta)
                    dF_lat = v_mag > 1e-20 ? dF * (v_lat / v_mag) : 0.0
                    dF_roll = v_mag > 1e-20 ? dF * (v_roll / v_mag) : 0.0

                    H_slide_o += dF_roll * u_slide_o_dim  # 【Bug8修复】有符号自然积分
                    H_spin_o += dF_lat * v_micro_lat + dF_roll * v_micro_roll
                end
            end
        end

        # ── 拖曳热 [W] (流体粘性 + Goksem-Hargreaves EHL 滚动阻力) ──
        H_drag = F_d_dim * V_ball_dim
        ν_cst = ρ_eff > 0 ? (μ_oil_dyn / ρ_eff) * 1e6 : 10.0
        ndm_factor = 60000.0 / π

        if Q_i_dim > 0 && abs(u_entrain_i) > 1e-6
            U_i = lub_mu0 * abs(u_entrain_i) / (E_prime_dim * 0.5 * D_ball_dim + 1e-30)
            G_i = lub_alpha_pv * E_prime_dim
            W_ehl_i = Q_i_dim / (E_prime_dim * (0.5 * D_ball_dim)^2 + 1e-30)
            if U_i * G_i > 0 && W_ehl_i > 0
                F_roll_i = 4.318 * E_prime_dim * (0.5 * D_ball_dim)^2 / G_i * (U_i * G_i)^0.658 * W_ehl_i^0.0126
                H_drag += F_roll_i * (1.0 / (1.0 + 1.84e-9 * (abs(u_entrain_i) * ndm_factor)^1.28 * ν_cst^0.64)) * abs(u_entrain_i)
            end
        end
        if Q_o_dim > 0 && abs(u_entrain_o) > 1e-6
            U_o = lub_mu0 * abs(u_entrain_o) / (E_prime_dim * 0.5 * D_ball_dim + 1e-30)
            G_o = lub_alpha_pv * E_prime_dim
            W_ehl_o = Q_o_dim / (E_prime_dim * (0.5 * D_ball_dim)^2 + 1e-30)
            if U_o * G_o > 0 && W_ehl_o > 0
                F_roll_o = 4.318 * E_prime_dim * (0.5 * D_ball_dim)^2 / G_o * (U_o * G_o)^0.658 * W_ehl_o^0.0126
                H_drag += F_roll_o * (1.0 / (1.0 + 1.84e-9 * (abs(u_entrain_o) * ndm_factor)^1.28 * ν_cst^0.64)) * abs(u_entrain_o)
            end
        end

        # ── 搅油热 [W] ──
        ω_ball_mag_dim = sqrt(ω_orb_x^2 + ω_orb_r^2 + ω_orb_θ^2) * W_scale
        M_ch_dim = churning_moment(ω_ball_mag_dim, 0.5 * D_star * L_scale, ρ_eff, μ_oil_dyn)
        H_churn = M_ch_dim * ω_ball_mag_dim

        # ── 写出场输出 (全部为有量纲标准单位) ──
        fo[base+FO_Q_I] = Q_i_dim                  # [N]
        fo[base+FO_Q_O] = Q_o_dim                  # [N]
        fo[base+FO_ALPHA_I] = α_i                   # [rad]
        fo[base+FO_ALPHA_O] = α_o                   # [rad]
        fo[base+FO_DELTA_I] = δ_i_sm * L_scale      # [m]
        fo[base+FO_DELTA_O] = δ_o_sm * L_scale      # [m]
        fo[base+FO_A_I] = a_i_dim                   # [m]
        fo[base+FO_B_I] = b_i_dim                   # [m]
        fo[base+FO_A_O] = a_o_dim                   # [m]
        fo[base+FO_B_O] = b_o_dim                   # [m]
        fo[base+FO_U_SLIDE_I] = u_slide_i_dim       # [m/s]
        fo[base+FO_U_SLIDE_O] = u_slide_o_dim       # [m/s]
        fo[base+FO_U_MEAN_I] = u_mean_i * V_scale   # [m/s]
        fo[base+FO_U_MEAN_O] = u_mean_o * V_scale   # [m/s]
        fo[base+FO_H_FILM_I] = h_film_i             # [m]
        fo[base+FO_H_FILM_O] = h_film_o             # [m]
        fo[base+FO_KAPPA_I] = κ_i                   # [-]
        fo[base+FO_KAPPA_O] = κ_o                   # [-]
        fo[base+FO_F_TRAC_I] = F_trac_i_dim         # [N]
        fo[base+FO_F_TRAC_O] = F_trac_o_dim         # [N]
        fo[base+FO_M_SPIN_I] = M_sp_i_dim           # [N·m]
        fo[base+FO_M_SPIN_O] = M_sp_o_dim           # [N·m]
        fo[base+FO_F_DRAG] = F_d_dim                 # [N]
        fo[base+FO_M_CHURN] = M_ch_dim               # [N·m]
        fo[base+FO_OMEGA_X] = ω_orb_x * W_scale     # [rad/s]
        fo[base+FO_OMEGA_Y] = ω_orb_r * W_scale     # [rad/s]
        fo[base+FO_OMEGA_Z] = ω_orb_θ * W_scale     # [rad/s]
        fo[base+FO_R_BALL] = r_b * L_scale           # [m]
        fo[base+FO_THETA_DOT] = θ̇_b * W_scale       # [rad/s]
        fo[base+FO_Q0] = real(q)
        fo[base+FO_Q1], fo[base+FO_Q2], fo[base+FO_Q3] = imag_part(q)
        fo[base+FO_H_SLIDE_I] = H_slide_i           # [W]
        fo[base+FO_H_SLIDE_O] = H_slide_o           # [W]
        fo[base+FO_H_SPIN_I] = H_spin_i             # [W]
        fo[base+FO_H_SPIN_O] = H_spin_o             # [W]
        fo[base+FO_H_DRAG] = H_drag                  # [W]
        fo[base+FO_H_CHURN] = H_churn               # [W]
        fo[base+FO_F_POCKET] = F_pk                  # [N]
        fo[base+FO_W_SPIN_I] = ω_spin_i_dim         # [rad/s]
        fo[base+FO_W_SPIN_O] = ω_spin_o_dim         # [rad/s]
    end

    return fo
end

# ── 辅助函数：Hamrock-Dowson EHL 膜厚 ──
"""
    _film_thickness_hd(u_mean, Q, E_prime, R_eff, kappa_e, mu_0, alpha_pv, K_th, T_0)

Hamrock-Dowson film thickness with thermal correction.
Returns film thickness h [m].
"""
function _film_thickness_hd(u_mean, Q, E_prime, R_eff, kappa_e, mu_0, alpha_pv, K_th, T_0)
    if Q <= 0 || abs(u_mean) < 1e-15 || R_eff <= 0 || E_prime <= 0
        return 0.0
    end
    U = mu_0 * abs(u_mean) / (E_prime * R_eff)
    G = alpha_pv * E_prime
    W = Q / (E_prime * R_eff^2)
    (U <= 0 || G <= 0 || W <= 0) && return 0.0

    ellip_factor = max(1.0 - 0.61 * exp(-0.73 * max(kappa_e, 0.1)), 0.01)

    ln_h_R = log(2.69) + 0.67 * log(U) + 0.53 * log(G) + log(ellip_factor) - 0.067 * log(W)
    h_iso = R_eff * exp(ln_h_R)

    phi_T = 1.0
    if K_th > 0 && T_0 > 0
        Q_m = 2.0 * mu_0 * u_mean^2 / (K_th * T_0)
        phi_T = Q_m > 1e-10 ? max(1.0 / (1.0 + 0.1 * Q_m^0.64), 0.01) : 1.0
    end

    return max(h_iso * phi_T, 1e-12)
end
