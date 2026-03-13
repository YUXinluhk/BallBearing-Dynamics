# =====================================================================
# Dynamics/kernel.jl — ODE RHS for 6-DOF bearing dynamics (究极完全体)
#
# 已补全所有强坐标耦合、兜孔三维摩擦、内圈陀螺力矩及反作用力配平
# =====================================================================

"""
    ode_rhs!(du, u, params_tuple, t)
"""
function ode_rhs!(du::AbstractVector, u::AbstractVector, params_tuple, t)

    p = params_tuple[1]
    Z = Int(p[P_NBALL])
    D_star = p[P_D]
    d_m_star = p[P_DM]
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
    β_temp = p[P_BETA_TEMP]
    Λ_LSS = p[P_LAMBDA_LSS]
    trac_A = p[P_TRAC_A]
    trac_B = p[P_TRAC_B]
    trac_C = p[P_TRAC_C]
    trac_D = p[P_TRAC_D]
    μ_spin = p[P_MU_SPIN]
    ζ = p[P_ZETA]
    F_a_star = p[P_FA]
    F_r_star = p[P_FR]
    t_ramp = p[P_T_RAMP]
    ρ_eff = p[P_RHO_EFF]
    μ_oil = p[P_MU_OIL]
    cage_web = p[P_CAGE_WEB]
    R_ball_k_nom = D_star / 2.0
    drb_i = p[P_DRB_I]
    drb_o = p[P_DRB_O]
    x_gi0_star = p[P_X_GI0]
    x_go0_star = p[P_X_GO0]
    f_i_geom = p[P_FI]
    f_o_geom = p[P_FO]

    # Damping & Scales
    c_ball_t = p[P_C_BALL_TRANS]
    c_ball_o = p[P_C_BALL_ORBIT]
    c_ball_spin = p[P_C_BALL_SPIN]
    c_ir_damp = p[P_C_IR_DAMP]
    c_cage_damp = p[P_C_CAGE_DAMP]
    c_tilt_damp = p[P_C_TILT]
    ω_cage_star = p[P_OMEGA_CAGE]
    V_scale = p[P_V_SCALE]
    W_scale = p[P_W_SCALE]
    L_scale = p[P_L_SCALE]
    Q_scale = p[P_Q_SCALE]
    d_dim = D_star * L_scale
    d_m_dim = d_m_star * L_scale
    E_prime_dim = p[P_E_PRIME] * Q_scale / L_scale^2
    power_scale = Q_scale * V_scale

    # Cage params
    cage_mass = p[P_CAGE_MASS]
    cage_Ixx = p[P_CAGE_IXX]
    k_pocket = p[P_K_POCKET]
    μ_pocket = p[P_MU_POCKET]
    c_cage = p[P_C_CAGE]

    # ── Load & speed ──
    ramp = t_ramp > 0 ? 0.5 * (1.0 - cos(π * min(t / t_ramp, 1.0))) : 1.0
    ω_ir = p[P_OMEGA_IR] * ramp
    ω_or = p[P_OMEGA_OR] * ramp
    effusivity = p[P_EFFUSIVITY]

    # ── Read state views ──
    ir_p = ir_pos_view(u, Z)
    ir_v = ir_vel_view(u, Z)
    x_ir, y_ir, z_ir, γy_ir, γz_ir = ir_p[1], ir_p[2], ir_p[3], ir_p[4], ir_p[5]
    ẋ_ir, ẏ_ir, ż_ir, γ̇y_ir, γ̇z_ir = ir_v[1], ir_v[2], ir_v[3], ir_v[4], ir_v[5]

    cg_p = cage_pos_view(u, Z)
    cg_v = cage_vel_view(u, Z)
    x_cg, y_cg, z_cg, θ_cg = cg_p[1], cg_p[2], cg_p[3], cg_p[4]
    ẋ_cg, ẏ_cg, ż_cg, θ̇_cg = cg_v[1], cg_v[2], cg_v[3], cg_v[4]

    tv = thermal_view(u, Z)
    T_ir_node, T_or_node, T_ball_node, T_oil_node = tv[1], tv[2], tv[3], tv[4]

    # Phase 4b: Ball thermal expansion — closes internal clearance at high DN
    δr_ball = p[P_CTE] * (T_ball_node - p[P_T_REF]) * R_ball_k_nom
    R_ball_k = R_ball_k_nom + δr_ball

    # ── Dynamic Viscosity (bidirectional Roelands: T>T_ref → μ↓, T<T_ref → μ↑) ──
    T_ref_oil = p[P_T0]
    μ_oil_dyn = μ_oil * exp(clamp(-β_temp * (T_oil_node - T_ref_oil), -30.0, 30.0))

    # ── Scalar accumulators (消除堆内存分配，赋予全局平滑阻尼) ──
    T_u = eltype(u)
    H_ir_total = zero(T_u)
    H_or_total = zero(T_u)
    H_ball_total = zero(T_u)
    H_oil_total = zero(T_u)

    F_ir_x = F_a_star - c_ir_damp * ẋ_ir
    F_ir_y = zero(T_u) - c_ir_damp * ẏ_ir
    F_ir_z = F_r_star - c_ir_damp * ż_ir
    M_ir_y = zero(T_u) - c_tilt_damp * γ̇y_ir
    M_ir_z = zero(T_u) - c_tilt_damp * γ̇z_ir

    F_cg_x = zero(T_u) - c_cage_damp * ẋ_cg
    F_cg_y = zero(T_u) - c_cage * ẏ_cg
    F_cg_z = zero(T_u) - c_cage * ż_cg
    F_cg_θ = zero(T_u) - c_cage_damp * (d_m_star / 2.0)^2 * (θ̇_cg - ω_cage_star)

    # ── 引导环支撑与摩擦闭环 (修复牛顿第三定律反作用力缺失) ──
    ecc_y = y_cg - (p[P_PILOT_IR] > 0.5 ? y_ir : zero(T_u))
    ecc_z = z_cg - (p[P_PILOT_IR] > 0.5 ? z_ir : zero(T_u))
    ecc = sqrt(ecc_y^2 + ecc_z^2 + 1e-30)

    δ_pilot = max(ecc - p[P_PILOT_CLR], zero(T_u))
    F_pilot = p[P_K_PILOT] * δ_pilot
    ey, ez = ecc_y / ecc, ecc_z / ecc

    F_cg_y -= F_pilot * ey
    F_cg_z -= F_pilot * ez

    Δω_pilot = θ̇_cg - (p[P_PILOT_IR] > 0.5 ? ω_ir : zero(T_u))
    R_pilot = p[P_PILOT_IR] > 0.5 ? p[P_CAGE_IR] : p[P_CAGE_OR]
    F_fric_pilot = p[P_MU_PILOT] * F_pilot * tanh(Δω_pilot / 1e-3)

    F_cg_θ -= F_fric_pilot * R_pilot
    F_cg_y += F_fric_pilot * ez
    F_cg_z -= F_fric_pilot * ey

    # 精确返还给内圈反作用力矩阵！
    if p[P_PILOT_IR] > 0.5
        F_ir_y += F_pilot * ey - F_fric_pilot * ez
        F_ir_z += F_pilot * ez + F_fric_pilot * ey
    end

    @inbounds for i in 1:N_IR_POS
        du[ir_pos_offset()+i-1] = ir_v[i]
    end
    @inbounds for i in 1:N_CAGE_POS
        du[cage_pos_offset(Z)+i-1] = cg_v[i]
    end

    # ── 滚珠主循环 ──
    @inbounds @fastmath for j in 1:Z
        bp = ball_pos_view(u, j, Z)
        bv = ball_vel_view(u, j, Z)
        x_b, r_b, θ_b = bp[1], bp[2], bp[3]
        # (ball_quat removed — quaternion accessed via raw u[] indices below for AD compatibility)
        ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]

        # 修正：定义严格右手系 ω_body = [x轴, r径向, θ切向]
        ω_body = ball_omega(u, j, Z)
        sθ, cθ = sincos(θ_b)

        # 【核心修复1】：深度几何干涉，内圈倾覆角不仅带来轴向位移，还带来径向杠杆偏心
        # Phase 4: 热膨胀动态游隙修正 — 参考温度为初始运行温度 T0
        # δr = CTE·(T - T_init)·R★ → t=0时δr=0 (与QS初值一致)
        δr_ir = p[P_CTE] * (T_ir_node - p[P_T_REF]) * (D_i_star / 2.0)
        δr_or = p[P_CTE] * (T_or_node - p[P_T_REF]) * (D_o_star / 2.0)
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

        cos_α_i = clamp(abs(nfi_r), 1e-4, 0.9999)
        cos_α_o = clamp(abs(nfo_r), 1e-4, 0.9999)
        Y_i_dim_rt, a_i_star_rt, b_i_star_rt, Σρ_i_dim_rt, E2_i_rt =
            hertz_runtime_contact(cos_α_i, d_dim, d_m_dim, f_i_geom, E_prime_dim, true)
        Y_o_dim_rt, a_o_star_rt, b_o_star_rt, Σρ_o_dim_rt, E2_o_rt =
            hertz_runtime_contact(cos_α_o, d_dim, d_m_dim, f_o_geom, E_prime_dim, false)

        Y_i_eff = Y_i_dim_rt > 0.0 ? Y_i_dim_rt * L_scale^1.5 / Q_scale : Y_i
        Y_o_eff = Y_o_dim_rt > 0.0 ? Y_o_dim_rt * L_scale^1.5 / Q_scale : Y_o
        a_i_star_eff = Y_i_dim_rt > 0.0 ? a_i_star_rt : a_i_star
        b_i_star_eff = Y_i_dim_rt > 0.0 ? b_i_star_rt : b_i_star
        a_o_star_eff = Y_o_dim_rt > 0.0 ? a_o_star_rt : a_o_star
        b_o_star_eff = Y_o_dim_rt > 0.0 ? b_o_star_rt : b_o_star
        Σρ_i_eff = Σρ_i_dim_rt > 1e-20 ? Σρ_i_dim_rt * L_scale : Σρ_i
        Σρ_o_eff = Σρ_o_dim_rt > 1e-20 ? Σρ_o_dim_rt * L_scale : Σρ_o
        E2_i_eff = Y_i_dim_rt > 0.0 ? E2_i_rt : p[P_E2_I]
        E2_o_eff = Y_o_dim_rt > 0.0 ? E2_o_rt : p[P_E2_O]

        δ_i_sm = smooth_hertz_delta(L_i - drb_i)
        δ_o_sm = smooth_hertz_delta(L_o - drb_o)
        Q_i = Y_i_eff * δ_i_sm * sqrt(δ_i_sm)
        Q_o = Y_o_eff * δ_o_sm * sqrt(δ_o_sm)

        c_i = Σρ_i_eff > 1e-20 ? cbrt(3Q_i / (2 * p[P_E_PRIME] * Σρ_i_eff)) : zero(T_u)
        a_i_dim, b_i_dim = a_i_star_eff * c_i * L_scale, b_i_star_eff * c_i * L_scale
        c_o = Σρ_o_eff > 1e-20 ? cbrt(3Q_o / (2 * p[P_E_PRIME] * Σρ_o_eff)) : zero(T_u)
        a_o_dim, b_o_dim = a_o_star_eff * c_o * L_scale, b_o_star_eff * c_o * L_scale

        # 3D contact geometry: predict α₂ from the lateral drift accumulated during
        # one contact traversal, then rebuild the contact basis with T_ac(α₁, α₂).
        e_x = SVector{3,T_u}(one(T_u), zero(T_u), zero(T_u))
        e_r = SVector{3,T_u}(zero(T_u), -sθ, cθ)
        e_θ = SVector{3,T_u}(zero(T_u), -cθ, -sθ)
        v_ball_θ = r_b * θ̇_b
        wx = ω_body[1]
        wy = -ω_body[2] * sθ - ω_body[3] * cθ
        wz = ω_body[2] * cθ - ω_body[3] * sθ
        ω_ball_global = SVector{3,T_u}(wx, wy, wz)
        v_ball_center_global = ẋ_b * e_x + ṙ_b * e_r + v_ball_θ * e_θ
        v_ir_center_global = SVector{3,T_u}(ẋ_ir, ẏ_ir, ż_ir)
        ω_ir_global = SVector{3,T_u}(ω_ir, γ̇y_ir, γ̇z_ir)
        ω_or_global = SVector{3,T_u}(ω_or, zero(T_u), zero(T_u))

        α1_i_0, _ = contact_angles_from_direction(dx_i, dr_i)
        α1_o_0, _ = contact_angles_from_direction(dx_o, dr_o)
        n_i_ac_0, t_lat_i_ac_0, _ = contact_basis_from_angles(α1_i_0, 0.0)
        n_o_ac_0, t_lat_o_ac_0, _ = contact_basis_from_angles(α1_o_0, 0.0)

        n_i_global_0 = n_i_ac_0[1] * e_x + n_i_ac_0[2] * e_θ + n_i_ac_0[3] * e_r
        n_o_global_0 = n_o_ac_0[1] * e_x + n_o_ac_0[2] * e_θ + n_o_ac_0[3] * e_r
        t_lat_i_global_0 = t_lat_i_ac_0[1] * e_x + t_lat_i_ac_0[2] * e_θ + t_lat_i_ac_0[3] * e_r
        t_lat_o_global_0 = t_lat_o_ac_0[1] * e_x + t_lat_o_ac_0[2] * e_θ + t_lat_o_ac_0[3] * e_r

        cp_i_global_0 = inertial_from_cylindrical(x_b, r_b, θ_b) - R_ball_k * n_i_global_0
        cp_o_global_0 = inertial_from_cylindrical(x_b, r_b, θ_b) - R_ball_k * n_o_global_0
        v_ball_cp_i_global_0 = v_ball_center_global + cross(ω_ball_global, -R_ball_k * n_i_global_0)
        v_ball_cp_o_global_0 = v_ball_center_global + cross(ω_ball_global, -R_ball_k * n_o_global_0)
        v_ir_cp_global_0 = v_ir_center_global + cross(ω_ir_global, cp_i_global_0 - SVector{3,T_u}(x_ir, y_ir, z_ir))
        v_or_cp_global_0 = cross(ω_or_global, cp_o_global_0)
        rel_i_global_0 = v_ir_cp_global_0 - v_ball_cp_i_global_0
        rel_o_global_0 = v_or_cp_global_0 - v_ball_cp_o_global_0
        u_side_i_pred = dot(rel_i_global_0, t_lat_i_global_0) * V_scale
        u_side_o_pred = dot(rel_o_global_0, t_lat_o_global_0) * V_scale

        v_pass_i = max(abs(dot(v_ir_cp_global_0 - v_ball_cp_i_global_0, e_θ)) * V_scale, 1e-6)
        v_pass_o = max(abs(dot(v_or_cp_global_0 - v_ball_cp_o_global_0, e_θ)) * V_scale, 1e-6)
        t_pass_i = max(a_i_dim, b_i_dim, 1e-9) / v_pass_i
        t_pass_o = max(a_o_dim, b_o_dim, 1e-9) / v_pass_o
        dθ_i_eff = clamp(u_side_i_pred * t_pass_i / L_scale, -0.25 * max(a_i_dim, b_i_dim) / L_scale, 0.25 * max(a_i_dim, b_i_dim) / L_scale)
        dθ_o_eff = clamp(u_side_o_pred * t_pass_o / L_scale, -0.25 * max(a_o_dim, b_o_dim) / L_scale, 0.25 * max(a_o_dim, b_o_dim) / L_scale)

        α1_i, α2_i = contact_angles_from_direction(dx_i, dθ_i_eff, dr_i)
        α1_o, α2_o = contact_angles_from_direction(dx_o, dθ_o_eff, dr_o)
        n_i_ac, t_lat_i_ac, t_roll_i_ac = contact_basis_from_angles(α1_i, α2_i)
        n_o_ac, t_lat_o_ac, t_roll_o_ac = contact_basis_from_angles(α1_o, α2_o)

        nfi_x, nfi_θ, nfi_r = n_i_ac[1], n_i_ac[2], n_i_ac[3]
        nfo_x, nfo_θ, nfo_r = n_o_ac[1], n_o_ac[2], n_o_ac[3]
        n_i_global = nfi_x * e_x + nfi_θ * e_θ + nfi_r * e_r
        n_o_global = nfo_x * e_x + nfo_θ * e_θ + nfo_r * e_r
        t_lat_i_x, t_lat_i_θ, t_lat_i_r = t_lat_i_ac[1], t_lat_i_ac[2], t_lat_i_ac[3]
        t_lat_o_x, t_lat_o_θ, t_lat_o_r = t_lat_o_ac[1], t_lat_o_ac[2], t_lat_o_ac[3]
        t_roll_i_global = t_roll_i_ac[1] * e_x + t_roll_i_ac[2] * e_θ + t_roll_i_ac[3] * e_r
        t_roll_o_global = t_roll_o_ac[1] * e_x + t_roll_o_ac[2] * e_θ + t_roll_o_ac[3] * e_r
        t_lat_i_global = t_lat_i_x * e_x + t_lat_i_θ * e_θ + t_lat_i_r * e_r
        t_lat_o_global = t_lat_o_x * e_x + t_lat_o_θ * e_θ + t_lat_o_r * e_r

        cos_α_i = clamp(abs(nfi_r), 1e-4, 0.9999)
        cos_α_o = clamp(abs(nfo_r), 1e-4, 0.9999)
        Y_i_dim_rt, a_i_star_rt, b_i_star_rt, Σρ_i_dim_rt, E2_i_rt =
            hertz_runtime_contact(cos_α_i, d_dim, d_m_dim, f_i_geom, E_prime_dim, true)
        Y_o_dim_rt, a_o_star_rt, b_o_star_rt, Σρ_o_dim_rt, E2_o_rt =
            hertz_runtime_contact(cos_α_o, d_dim, d_m_dim, f_o_geom, E_prime_dim, false)

        Y_i_eff = Y_i_dim_rt > 0.0 ? Y_i_dim_rt * L_scale^1.5 / Q_scale : Y_i
        Y_o_eff = Y_o_dim_rt > 0.0 ? Y_o_dim_rt * L_scale^1.5 / Q_scale : Y_o
        a_i_star_eff = Y_i_dim_rt > 0.0 ? a_i_star_rt : a_i_star
        b_i_star_eff = Y_i_dim_rt > 0.0 ? b_i_star_rt : b_i_star
        a_o_star_eff = Y_o_dim_rt > 0.0 ? a_o_star_rt : a_o_star
        b_o_star_eff = Y_o_dim_rt > 0.0 ? b_o_star_rt : b_o_star
        Σρ_i_eff = Σρ_i_dim_rt > 1e-20 ? Σρ_i_dim_rt * L_scale : Σρ_i
        Σρ_o_eff = Σρ_o_dim_rt > 1e-20 ? Σρ_o_dim_rt * L_scale : Σρ_o
        E2_i_eff = Y_i_dim_rt > 0.0 ? E2_i_rt : p[P_E2_I]
        E2_o_eff = Y_o_dim_rt > 0.0 ? E2_o_rt : p[P_E2_O]

        δ_i_sm = smooth_hertz_delta(sqrt(dx_i^2 + dθ_i_eff^2 + dr_i^2 + 1e-30) - drb_i)
        δ_o_sm = smooth_hertz_delta(sqrt(dx_o^2 + dθ_o_eff^2 + dr_o^2 + 1e-30) - drb_o)
        Q_i = Y_i_eff * δ_i_sm * sqrt(δ_i_sm)
        Q_o = Y_o_eff * δ_o_sm * sqrt(δ_o_sm)

        c_i = Σρ_i_eff > 1e-20 ? cbrt(3Q_i / (2 * p[P_E_PRIME] * Σρ_i_eff)) : zero(T_u)
        a_i_dim, b_i_dim = a_i_star_eff * c_i * L_scale, b_i_star_eff * c_i * L_scale
        c_o = Σρ_o_eff > 1e-20 ? cbrt(3Q_o / (2 * p[P_E_PRIME] * Σρ_o_eff)) : zero(T_u)
        a_o_dim, b_o_dim = a_o_star_eff * c_o * L_scale, b_o_star_eff * c_o * L_scale

        cp_i_global = inertial_from_cylindrical(x_b, r_b, θ_b) - R_ball_k * n_i_global
        cp_o_global = inertial_from_cylindrical(x_b, r_b, θ_b) - R_ball_k * n_o_global
        v_ball_cp_i_global = v_ball_center_global + cross(ω_ball_global, -R_ball_k * n_i_global)
        v_ball_cp_o_global = v_ball_center_global + cross(ω_ball_global, -R_ball_k * n_o_global)
        v_ir_cp_global = v_ir_center_global + cross(ω_ir_global, cp_i_global - SVector{3,T_u}(x_ir, y_ir, z_ir))
        v_or_cp_global = cross(ω_or_global, cp_o_global)
        rel_i_global = v_ir_cp_global - v_ball_cp_i_global
        rel_o_global = v_or_cp_global - v_ball_cp_o_global
        x_contact_i = cp_i_global[1]
        x_contact_o = cp_o_global[1]
        r_contact_i = dot(cp_i_global, e_r)
        r_contact_o = dot(cp_o_global, e_r)
        v_ir_roll = dot(v_ir_cp_global, t_roll_i_global)
        v_or_roll = dot(v_or_cp_global, t_roll_o_global)
        v_ball_surface_i = dot(v_ball_cp_i_global, t_roll_i_global)
        v_ball_surface_o = dot(v_ball_cp_o_global, t_roll_o_global)

        u_slide_i_dim = dot(rel_i_global, t_roll_i_global) * V_scale
        u_slide_o_dim = dot(rel_o_global, t_roll_o_global) * V_scale
        u_side_i_dim = dot(rel_i_global, t_lat_i_global) * V_scale
        u_side_o_dim = dot(rel_o_global, t_lat_o_global) * V_scale

        ω_ir_dim = ω_ir * W_scale
        ω_ir_spin_dim = ω_ir_dim * nfi_x + (-γ̇y_ir * sθ + γ̇z_ir * cθ) * W_scale * nfi_r
        ω_spin_i = ω_ir_spin_dim - (ω_body[1] * nfi_x + ω_body[2] * nfi_r) * W_scale
        ω_or_spin_dim = ω_or * W_scale * nfo_x
        ω_spin_o = ω_or_spin_dim - (ω_body[1] * nfo_x + ω_body[2] * nfo_r) * W_scale

        # 引入微观热弹流 (TEHD) 的抗发散牵引模型
        P_mean_i = Q_i > 0 ? (Q_i * Q_scale) / (π * a_i_dim * b_i_dim + 1e-16) : zero(T_u)
        P_mean_o = Q_o > 0 ? (Q_o * Q_scale) / (π * a_o_dim * b_o_dim + 1e-16) : zero(T_u)
        # 【Bug1修复】EHL卷吸速度：必须在接触点共转参考系下计算
        v_cp_i_dim = r_contact_i * θ̇_b * V_scale  # 内圈接触点轨道速度
        v_cp_o_dim = r_contact_o * θ̇_b * V_scale  # 外圈接触点轨道速度
        u_entrain_i = abs(0.5 * (v_ir_θ * V_scale + v_ball_surface_i * V_scale) - v_cp_i_dim)
        u_entrain_o = abs(0.5 * (v_or_θ * V_scale + v_ball_surface_o * V_scale) - v_cp_o_dim)

        ω_roll_abs_dim = sqrt(ω_body[1]^2 + ω_body[2]^2 + ω_body[3]^2 + 1e-16) * W_scale
        F_trac_i_dim, F_side_i_dim, H_slide_i_dim, H_spin_i_dim = integrate_tehd_contact_force(
            u_slide_i_dim, u_side_i_dim, P_mean_i, a_i_dim, b_i_dim, u_entrain_i,
            ω_spin_i, ω_roll_abs_dim, R_ball_k * L_scale, f_i_geom,
            trac_A, trac_B, trac_C, trac_D, Λ_LSS, β_temp, effusivity;
            flash_length=b_i_dim,
        )
        F_trac_o_dim, F_side_o_dim, H_slide_o_dim, H_spin_o_dim = integrate_tehd_contact_force(
            u_slide_o_dim, u_side_o_dim, P_mean_o, a_o_dim, b_o_dim, u_entrain_o,
            ω_spin_o, ω_roll_abs_dim, R_ball_k * L_scale, f_o_geom,
            trac_A, trac_B, trac_C, trac_D, Λ_LSS, β_temp, effusivity;
            flash_length=0.5 * a_o_dim,
        )
        F_trac_tang_i = F_trac_i_dim / Q_scale
        F_trac_tang_o = F_trac_o_dim / Q_scale
        F_trac_lat_i = F_side_i_dim / Q_scale
        F_trac_lat_o = F_side_o_dim / Q_scale

        # 自旋摩擦力矩 (Harris 完整公式, 含椭圆积分 E(m))
        M_spin_i = spin_moment(ω_spin_i, Q_i, a_i_star_eff * c_i, b_i_star_eff * c_i, μ_spin, E2_i_eff)
        M_spin_o = spin_moment(ω_spin_o, Q_o, a_o_star_eff * c_o, b_o_star_eff * c_o, μ_spin, E2_o_eff)

        # ── 强耦合热网络发热累加 (Non-Dimensional Power) ──
        P_slide_i_nd = H_slide_i_dim / power_scale
        P_slide_o_nd = H_slide_o_dim / power_scale
        P_spin_i_nd = H_spin_i_dim / power_scale
        P_spin_o_nd = H_spin_o_dim / power_scale

        H_i_tot = P_slide_i_nd + P_spin_i_nd
        H_o_tot = P_slide_o_nd + P_spin_o_nd

        H_ir_total += 0.5 * H_i_tot
        H_or_total += 0.5 * H_o_tot
        H_ball_total += 0.5 * (H_i_tot + H_o_tot)

        # ── 3D 流体力学 ──
        v_rel_θ_nd = r_b * (θ̇_b - ω_cage_star)
        V_ball_dim = sqrt(ẋ_b^2 + v_rel_θ_nd^2 + ṙ_b^2 + 1e-30) * V_scale
        F_drag_total = drag_force(V_ball_dim, D_star * L_scale, ρ_eff, μ_oil_dyn, cage_web * L_scale) / Q_scale
        F_drag_x = -F_drag_total * (ẋ_b * V_scale / V_ball_dim)
        F_drag_r = -F_drag_total * (ṙ_b * V_scale / V_ball_dim)
        F_drag_θ = -F_drag_total * (v_rel_θ_nd * V_scale / V_ball_dim)

        ω_mag_dim = sqrt(ω_body[1]^2 + ω_body[2]^2 + ω_body[3]^2 + 1e-30) * W_scale
        M_ch_nondim = churning_moment(ω_mag_dim, R_ball_k * L_scale, ρ_eff, μ_oil_dyn) / (Q_scale * L_scale)

        P_drag_nd = F_drag_total * (V_ball_dim / V_scale)
        P_drag_spin_nd = M_ch_nondim * (ω_mag_dim / W_scale)
        H_oil_total += P_drag_nd + P_drag_spin_nd

        # 【核心修复2】：保持架兜孔推力完美分解入三维坐标系，触发真实偏心涡动！
        pocket_θ = θ_cg + (j - 1) * 2π / Z
        sψ, cψ = sincos(pocket_θ)
        Δθ_pk = θ_b - pocket_θ
        Δθ_pk -= 2π * round(Δθ_pk / (2π))

        pen_pk_sm = smooth_hertz_delta(abs(r_b * Δθ_pk) - p[P_POCKET_CLR])
        F_pk_norm = p[P_K_POCKET] * pen_pk_sm

        F_pc_tang_ball = -F_pk_norm * tanh((r_b * Δθ_pk) * L_scale / 1e-5)

        # 兜孔内干摩擦高频耗散，驱动保持架轴向窜动
        v_cg_r_local = -ẏ_cg * sψ + ż_cg * cψ
        # 【Bug6修复】2D库仑摩擦：先合成滑移矢量，再投影方向，避免正方形极限域
        v_slip_ax = (ẋ_b - ẋ_cg) * V_scale
        v_slip_rad = (ṙ_b - v_cg_r_local) * V_scale
        v_slip_mag = sqrt(v_slip_ax^2 + v_slip_rad^2 + 1e-12)
        fric_mag = μ_pocket * F_pk_norm * tanh(v_slip_mag / 0.01)
        F_pc_ax_ball = -fric_mag * (v_slip_ax / v_slip_mag)
        F_pc_rad_ball = -fric_mag * (v_slip_rad / v_slip_mag)

        F_ball_x = Q_i * nfi_x + Q_o * nfo_x +
                   F_trac_lat_i * t_lat_i_x + F_trac_lat_o * t_lat_o_x +
                   F_drag_x + F_pc_ax_ball - c_ball_t * ẋ_b
        F_ball_r = Q_i * nfi_r + Q_o * nfo_r +
                   F_trac_lat_i * t_lat_i_r + F_trac_lat_o * t_lat_o_r +
                   F_drag_r + F_pc_rad_ball - c_ball_t * ṙ_b
        F_ball_θ = F_trac_tang_i + F_trac_tang_o + F_drag_θ + F_pc_tang_ball - c_ball_o * r_b * (θ̇_b - ω_cage_star)

        # 将兜孔推力完美映射到全局笛卡尔保持架方程！
        F_cg_θ -= F_pc_tang_ball * r_b
        F_cg_y += F_pc_rad_ball * sψ + F_pc_tang_ball * cψ
        F_cg_z += -F_pc_rad_ball * cψ + F_pc_tang_ball * sψ
        F_cg_x -= F_pc_ax_ball

        # ── 3D 扭矩体系与右手系四元数积分 ──
        moment_scale_inv = 1.0 / (Q_scale * L_scale)
        R_ball_dim = R_ball_k * L_scale
        τ_roll_i_x = -R_ball_dim * nfi_r * (F_trac_tang_i * Q_scale) * moment_scale_inv
        τ_roll_i_r = R_ball_dim * nfi_x * (F_trac_tang_i * Q_scale) * moment_scale_inv
        τ_roll_o_x = -R_ball_dim * nfo_r * (F_trac_tang_o * Q_scale) * moment_scale_inv
        τ_roll_o_r = R_ball_dim * nfo_x * (F_trac_tang_o * Q_scale) * moment_scale_inv

        # 【核心修复9】3D 兜孔摩擦力矩完美闭环
        # 【Bug5修复】3D 兜孔摩擦力矩：严格右手系叉乘 r×F
        τ_pk_x = sign(Δθ_pk) * R_ball_k * F_pc_rad_ball
        τ_pk_r = -sign(Δθ_pk) * R_ball_k * F_pc_ax_ball

        M_ball_x = τ_roll_i_x + τ_roll_o_x + M_spin_i * nfi_x + M_spin_o * nfo_x + τ_pk_x - c_ball_spin * ω_body[1] - M_ch_nondim * (ω_body[1] * W_scale / ω_mag_dim)
        M_ball_r = τ_roll_i_r + τ_roll_o_r + M_spin_i * nfi_r + M_spin_o * nfo_r + τ_pk_r - c_ball_spin * ω_body[2] - M_ch_nondim * (ω_body[2] * W_scale / ω_mag_dim)
        M_ball_θ = -R_ball_k * (F_trac_lat_i + F_trac_lat_o) - c_ball_spin * ω_body[3] - M_ch_nondim * (ω_body[3] * W_scale / ω_mag_dim)

        # 【核心修复3】：欧拉进动角动量修正 (-Ω × H) 恢复自然进动！
        M_ball_r += J_ball * θ̇_b * ω_body[3]
        M_ball_θ -= J_ball * θ̇_b * ω_body[2]

        off_p = ball_pos_offset(j)
        du[off_p], du[off_p+1], du[off_p+2] = ẋ_b, ṙ_b, θ̇_b

        # 【核心修复4】：原生的泛型化四元数逆变换与积分，解除包约束，赋能 ForwardDiff！
        wx, wy, wz = ω_body[1], -ω_body[2] * sθ - ω_body[3] * cθ, ω_body[2] * cθ - ω_body[3] * sθ
        qw, qx, qy, qz = u[off_p+3], u[off_p+4], u[off_p+5], u[off_p+6]

        # inv_rotate_vector 展开与积分，已移至 Transfers/quaternion.jl 中以解耦和复用
        du[off_p+3], du[off_p+4], du[off_p+5], du[off_p+6] = kinematics_quat_derivative_baumgarte(qw, qx, qy, qz, wx, wy, wz, 50.0)

        off_v = ball_vel_offset(j, Z)
        du[off_v], du[off_v+1], du[off_v+2] = F_ball_x / m_ball, F_ball_r / m_ball + r_b * θ̇_b^2, (F_ball_θ / m_ball - 2.0 * ṙ_b * θ̇_b) / (r_b + 1e-30)
        du[off_v+3], du[off_v+4], du[off_v+5] = M_ball_x / J_ball, M_ball_r / J_ball, M_ball_θ / J_ball

        # 【核心修复5】：采用精确物理触点作为力臂，构建毫无遗漏的内圈反作用张量矩阵！
        F_ir_X_add = -Q_i * nfi_x - F_trac_lat_i * t_lat_i_x
        F_ir_Y_add = Q_i * nfi_r * sθ + F_trac_tang_i * cθ - F_trac_lat_i * t_lat_i_r * (-sθ)
        F_ir_Z_add = -Q_i * nfi_r * cθ + F_trac_tang_i * sθ - F_trac_lat_i * t_lat_i_r * cθ

        F_ir_x += F_ir_X_add
        F_ir_y += F_ir_Y_add
        F_ir_z += F_ir_Z_add

        X_rel = x_contact_i - x_ir
        Y_rel = -r_contact_i * sθ - y_ir
        Z_rel = r_contact_i * cθ - z_ir

        M_ir_y += Z_rel * F_ir_X_add - X_rel * F_ir_Z_add + M_spin_i * nfi_r * sθ
        M_ir_z += X_rel * F_ir_Y_add - Y_rel * F_ir_X_add - M_spin_i * nfi_r * cθ
    end

    # ── 【核心修复6】内圈超强陀螺刚化力矩 (Gyroscopic Stiffening) 闭环 ──
    I_p_ir = m_ir * (D_i_star / 2)^2    # 极转动惯量
    I_d_ir = I_p_ir / 2.0               # 径向倾覆惯量

    off_irv = ir_vel_offset(Z)
    du[off_irv] = F_ir_x / m_ir
    du[off_irv+1] = F_ir_y / m_ir
    du[off_irv+2] = F_ir_z / m_ir
    du[off_irv+3] = (M_ir_y - I_p_ir * ω_ir * γ̇z_ir) / I_d_ir
    du[off_irv+4] = (M_ir_z + I_p_ir * ω_ir * γ̇y_ir) / I_d_ir

    off_cv = cage_vel_offset(Z)
    du[off_cv] = F_cg_x / cage_mass
    du[off_cv+1] = F_cg_y / cage_mass
    du[off_cv+2] = F_cg_z / cage_mass
    du[off_cv+3] = F_cg_θ / cage_Ixx

    # ── 【核心修复7】热网络强耦合：将无量纲热量转化为各节点温升 ──
    # Q_cond = G * (T1 - T2) in dimensionless form
    T_amb = p[P_T_AMB]

    Q_ir_ball = p[P_G_IR_BALL] * (T_ir_node - T_ball_node)
    Q_or_ball = p[P_G_OR_BALL] * (T_or_node - T_ball_node)
    Q_ball_oil = p[P_G_BALL_OIL] * (T_ball_node - T_oil_node)
    Q_or_amb = p[P_G_OR_AMB] * (T_or_node - T_amb)
    Q_oil_amb = p[P_G_OIL_AMB] * (T_oil_node - T_amb)

    # Oil flow circulation cooling: Q_flow = ṁ·cₚ·(T_oil - T_oil_inlet)
    # When oil_flow_rate > 0, this represents convective heat removal by flowing lubricant
    Q_oil_flow = p[P_OIL_FLOW_MDOT_CP] * (T_oil_node - p[P_T_OIL_INLET])

    # Note: T_ball_node is the single lumped temperature of ALL balls, and H_ball_total is the heat generated by ALL balls.
    dT_ir_dt = (H_ir_total - Q_ir_ball) / p[P_MCP_I]
    dT_or_dt = (H_or_total - Q_or_ball - Q_or_amb) / p[P_MCP_O]
    dT_ball_dt = (H_ball_total + Q_ir_ball + Q_or_ball - Q_ball_oil) / p[P_MCP_BALL]
    dT_oil_dt = (H_oil_total + Q_ball_oil - Q_oil_amb - Q_oil_flow) / p[P_MCP_OIL]

    off_th = thermal_offset(Z)
    du[off_th] = dT_ir_dt
    du[off_th+1] = dT_or_dt
    du[off_th+2] = dT_ball_dt
    du[off_th+3] = dT_oil_dt

    # ── 热量积分器: dE/dt = H → ODE 求解器精确积分 ──
    off_ha = heat_accum_offset(Z)
    du[off_ha] = H_ir_total                  # ∫H_ir dτ
    du[off_ha+1] = H_or_total                  # ∫H_or dτ
    du[off_ha+2] = H_ball_total                # ∫H_ball dτ
    du[off_ha+3] = H_oil_total                 # ∫H_oil dτ
    du[off_ha+4] = Q_or_amb + Q_oil_amb        # ∫Q_to_amb dτ
    # ∫Σ C★·dT/dτ dτ: ODE-consistent stored energy (avoids FBDF Newton residual drift)
    du[off_ha+5] = p[P_MCP_I] * dT_ir_dt + p[P_MCP_O] * dT_or_dt + p[P_MCP_BALL] * dT_ball_dt + p[P_MCP_OIL] * dT_oil_dt

    return nothing
end
