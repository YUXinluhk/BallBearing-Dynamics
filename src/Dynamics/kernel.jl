# =====================================================================
# Dynamics/kernel.jl — ODE RHS for 6-DOF bearing dynamics (究极完全体)
#
# 已补全所有强坐标耦合、兜孔三维摩擦、内圈陀螺力矩及反作用力配平
# =====================================================================

"""
    ode_rhs!(du, u, params_tuple, t)
"""
function ode_rhs!(du::AbstractVector, u::AbstractVector, params_tuple, t)

    p = params_tuple[1]::ODEParams
    Z = p.geom.Z
    D_star = p.geom.D
    d_m_star = p.geom.d_m
    m_ball = p.mass.m_ball
    J_ball = p.mass.J_ball
    m_ir = p.mass.m_ir
    Y_i = p.hertz.Y_i
    Y_o = p.hertz.Y_o
    a_i_star = p.hertz.a_i
    b_i_star = p.hertz.b_i
    a_o_star = p.hertz.a_o
    b_o_star = p.hertz.b_o
    Σρ_i = p.hertz.sr_i
    Σρ_o = p.hertz.sr_o
    D_i_star = p.geom.D_i
    D_o_star = p.geom.D_o
    μ₀ = p.lub.mu_0
    α_pv = p.lub.alpha_pv
    β_temp = p.lub.beta_temp
    Λ_LSS = p.lub.Lambda_LSS
    trac_A = p.trac.A
    trac_B = p.trac.B
    trac_C = p.trac.C
    trac_D = p.trac.D
    μ_spin = p.load.mu_spin
    ζ = p.load.zeta
    F_a_star = p.load.F_a
    F_r_star = p.load.F_r
    t_ramp = p.load.t_ramp
    ρ_eff = p.lub.rho_eff
    μ_oil = p.lub.mu_oil
    cage_web = p.cage.cage_web
    R_ball_k_nom = D_star / 2.0
    drb_i = p.geom.drb_i
    drb_o = p.geom.drb_o
    x_gi0_star = p.geom.x_gi0
    x_go0_star = p.geom.x_go0
    f_i_geom = p.geom.f_i
    f_o_geom = p.geom.f_o

    # Damping & Scales
    c_ball_t = p.damp.c_ball_trans
    c_ball_o = p.damp.c_ball_orbit
    c_ball_spin = p.damp.c_ball_spin
    c_ir_damp = p.damp.c_ir
    c_cage_damp = p.damp.c_cage
    c_tilt_damp = p.damp.c_tilt
    ω_cage_star = p.load.omega_cage
    V_scale = p.scale.V
    W_scale = p.scale.W
    L_scale = p.scale.L
    Q_scale = p.scale.Q
    d_dim = D_star * L_scale
    d_m_dim = d_m_star * L_scale
    E_prime_dim = p.hertz.E_prime * Q_scale / L_scale^2

    # Cage params
    cage_mass = p.mass.cage_mass
    cage_Ixx = p.mass.cage_Ixx
    k_pocket = p.cage.k_pocket
    μ_pocket = p.cage.mu_pocket
    c_cage = p.cage.c_cage

    ctx = build_kinematic_context(u, p, t)
    
    T_u = eltype(u)
    H_ir_total = zero(T_u)
    H_or_total = zero(T_u)
    H_ball_total = zero(T_u)
    H_oil_total = zero(T_u)

    F_ir_x = F_a_star - c_ir_damp * ctx.ẋ_ir
    F_ir_y = zero(T_u) - c_ir_damp * ctx.ẏ_ir
    F_ir_z = F_r_star - c_ir_damp * ctx.ż_ir
    M_ir_y = zero(T_u) - c_tilt_damp * ctx.γ̇y_ir
    M_ir_z = zero(T_u) - c_tilt_damp * ctx.γ̇z_ir

    F_cg_x = zero(T_u) - c_cage_damp * ctx.ẋ_cg
    F_cg_y = zero(T_u) - c_cage * ctx.ẏ_cg
    F_cg_z = zero(T_u) - c_cage * ctx.ż_cg
    F_cg_θ = zero(T_u) - c_cage_damp * (d_m_star / 2.0)^2 * (ctx.θ̇_cg - ω_cage_star)

    # ── 引导环支撑与摩擦闭环 ──
    ecc_y = ctx.y_cg - (p.cage.pilot_is_inner ? ctx.y_ir : zero(T_u))
    ecc_z = ctx.z_cg - (p.cage.pilot_is_inner ? ctx.z_ir : zero(T_u))
    ecc = sqrt(ecc_y^2 + ecc_z^2 + 1e-30)

    δ_pilot = max(ecc - p.cage.pilot_clr, zero(T_u))
    F_pilot = p.cage.k_pilot * δ_pilot
    ey, ez = ecc_y / ecc, ecc_z / ecc

    F_cg_y -= F_pilot * ey
    F_cg_z -= F_pilot * ez

    Δω_pilot = ctx.θ̇_cg - (p.cage.pilot_is_inner ? ctx.ω_ir : zero(T_u))
    R_pilot = p.cage.pilot_is_inner ? p.cage.cage_ir : p.cage.cage_or
    F_fric_pilot = p.cage.mu_pilot * F_pilot * tanh(Δω_pilot / 1e-3)

    F_cg_θ -= F_fric_pilot * R_pilot
    F_cg_y += F_fric_pilot * ez
    F_cg_z -= F_fric_pilot * ey

    if p.cage.pilot_is_inner
        F_ir_y += F_pilot * ey - F_fric_pilot * ez
        F_ir_z += F_pilot * ez + F_fric_pilot * ey
    end

    # ── IR & cage position derivatives = velocity (kinematic identity) ──
    du.ir.pos.x  = ctx.ẋ_ir
    du.ir.pos.y  = ctx.ẏ_ir
    du.ir.pos.z  = ctx.ż_ir
    du.ir.pos.γy = ctx.γ̇y_ir
    du.ir.pos.γz = ctx.γ̇z_ir

    du.cage.pos.x = ctx.ẋ_cg
    du.cage.pos.y = ctx.ẏ_cg
    du.cage.pos.z = ctx.ż_cg
    du.cage.pos.θ = ctx.θ̇_cg

    # ── 滚珠主循环 ──
    @inbounds for j in 1:Z
        res = compute_ball_tehd_kinematics(ctx, j)
        
        H_ir_total += 0.5 * res.H_i_tot
        H_or_total += 0.5 * res.H_o_tot
        H_ball_total += 0.5 * (res.H_i_tot + res.H_o_tot)
        H_oil_total += res.P_drag_nd + res.P_drag_spin_nd + res.P_ehl_roll_nd

        bp = u.ball[j].pos
        bv = u.ball[j].vel
        x_b, r_b, θ_b = bp.x, bp.r, bp.θ
        ẋ_b, ṙ_b, θ̇_b = bv.x, bv.r, bv.θ
        ω_body = ball_omega(u, j, Z)
        sθ, cθ = sincos(θ_b)

        pocket_θ = ctx.θ_cg + (j - 1) * 2π / Z
        sψ, cψ = sincos(pocket_θ)
        Δθ_pk = θ_b - pocket_θ
        Δθ_pk -= 2π * round(Δθ_pk / (2π))

        pen_pk_sm = smooth_hertz_delta(hypot(r_b * Δθ_pk, T_u(1e-12)) - p.cage.pocket_clr; ε=p.hertz.eps_contact)
        F_pk_norm = p.cage.k_pocket * pen_pk_sm
        F_pc_tang_ball = -F_pk_norm * tanh((r_b * Δθ_pk) * ctx.L_scale / 1e-5)

        v_cg_r_local = -ctx.ẏ_cg * sψ + ctx.ż_cg * cψ
        v_slip_ax = (ẋ_b - ctx.ẋ_cg) * ctx.V_scale
        v_slip_rad = (ṙ_b - v_cg_r_local) * ctx.V_scale
        v_slip_mag = sqrt(v_slip_ax^2 + v_slip_rad^2 + 1e-12)
        fric_mag = μ_pocket * F_pk_norm * tanh(v_slip_mag / 0.01)
        F_pc_ax_ball = -fric_mag * (v_slip_ax / v_slip_mag)
        F_pc_rad_ball = -fric_mag * (v_slip_rad / v_slip_mag)

        F_contact_i_loc = res.Q_i * res.n_i_loc + res.F_side_i * res.t_lat_i_loc + res.F_trac_tang_i * res.t_roll_i_loc
        F_contact_o_loc = res.Q_o * res.n_o_loc + res.F_side_o * res.t_lat_o_loc + res.F_trac_tang_o * res.t_roll_o_loc
        
        F_ball_x = F_contact_i_loc[1] + F_contact_o_loc[1] + res.F_drag_x + F_pc_ax_ball - c_ball_t * ẋ_b
        F_ball_r = F_contact_i_loc[2] + F_contact_o_loc[2] + res.F_drag_r + F_pc_rad_ball - c_ball_t * ṙ_b
        F_ehl_retard_θ = -(res.F_roll_ehl_i_nd + res.F_roll_ehl_o_nd) * tanh(θ̇_b / T_u(1e-3))
        F_ball_θ = F_contact_i_loc[3] + F_contact_o_loc[3] + res.F_drag_θ + F_pc_tang_ball + F_ehl_retard_θ - c_ball_o * r_b * (θ̇_b - ω_cage_star)

        F_cg_θ -= F_pc_tang_ball * r_b
        F_cg_y += F_pc_rad_ball * sψ + F_pc_tang_ball * cψ
        F_cg_z += -F_pc_rad_ball * cψ + F_pc_tang_ball * sψ
        F_cg_x -= F_pc_ax_ball

        smooth_sign_pk = tanh((r_b * Δθ_pk) * ctx.L_scale / 1e-5)
        τ_pk_x = -smooth_sign_pk * res.R_ball_k * F_pc_rad_ball
        τ_pk_r =  smooth_sign_pk * res.R_ball_k * F_pc_ax_ball
        
        τ_contact_i_loc = res.R_ball_k * (res.F_trac_tang_i * res.t_lat_i_loc - res.F_side_i * res.t_roll_i_loc)
        τ_contact_o_loc = res.R_ball_k * (res.F_trac_tang_o * res.t_lat_o_loc - res.F_side_o * res.t_roll_o_loc)
        
        M_ball_x = τ_contact_i_loc[1] + τ_contact_o_loc[1] + res.M_spin_i_nd * res.n_i_loc[1] + res.M_spin_o_nd * res.n_o_loc[1] + τ_pk_x - c_ball_spin * ω_body[1] - res.M_ch_nondim * (ω_body[1] * ctx.W_scale / sqrt(ω_body[1]^2 + ω_body[2]^2 + ω_body[3]^2 + 1e-30))
        M_ball_r = τ_contact_i_loc[2] + τ_contact_o_loc[2] + res.M_spin_i_nd * res.n_i_loc[2] + res.M_spin_o_nd * res.n_o_loc[2] + τ_pk_r - c_ball_spin * ω_body[2] - res.M_ch_nondim * (ω_body[2] * ctx.W_scale / sqrt(ω_body[1]^2 + ω_body[2]^2 + ω_body[3]^2 + 1e-30))
        M_ball_θ = τ_contact_i_loc[3] + τ_contact_o_loc[3] + res.M_spin_i_nd * res.n_i_loc[3] + res.M_spin_o_nd * res.n_o_loc[3] - c_ball_spin * ω_body[3] - res.M_ch_nondim * (ω_body[3] * ctx.W_scale / sqrt(ω_body[1]^2 + ω_body[2]^2 + ω_body[3]^2 + 1e-30))

        M_ball_r += J_ball * θ̇_b * ω_body[3]
        M_ball_θ -= J_ball * θ̇_b * ω_body[2]

        # ── Ball position derivatives (kinematic identity) ──
        du.ball[j].pos.x  = ẋ_b
        du.ball[j].pos.r  = ṙ_b
        du.ball[j].pos.θ  = θ̇_b

        wx, wy, wz = ω_body[1], -ω_body[2] * sθ - ω_body[3] * cθ, ω_body[2] * cθ - ω_body[3] * sθ
        bp_ca = u.ball[j].pos
        qw, qx, qy, qz = bp_ca.q0, bp_ca.q1, bp_ca.q2, bp_ca.q3

        dq0, dq1, dq2, dq3 = kinematics_quat_derivative_baumgarte(qw, qx, qy, qz, wx, wy, wz, 50.0)
        du.ball[j].pos.q0 = dq0
        du.ball[j].pos.q1 = dq1
        du.ball[j].pos.q2 = dq2
        du.ball[j].pos.q3 = dq3

        # ── Ball velocity derivatives (Newton-Euler in rotating frame) ──
        du.ball[j].vel.x  = F_ball_x / m_ball
        du.ball[j].vel.r  = F_ball_r / m_ball + r_b * θ̇_b^2
        du.ball[j].vel.θ  = (F_ball_θ / m_ball - 2.0 * ṙ_b * θ̇_b) / (r_b + 1e-30)
        du.ball[j].vel.ωx = M_ball_x / J_ball
        du.ball[j].vel.ωy = M_ball_r / J_ball
        du.ball[j].vel.ωz = M_ball_θ / J_ball

        F_ir_loc = -F_contact_i_loc
        F_ir_X_add = F_ir_loc[1]
        F_ir_Y_add = -F_ir_loc[2] * sθ - F_ir_loc[3] * cθ
        F_ir_Z_add =  F_ir_loc[2] * cθ - F_ir_loc[3] * sθ

        F_ir_x += F_ir_X_add
        F_ir_y += F_ir_Y_add
        F_ir_z += F_ir_Z_add

        M_ir_spin_loc = -res.M_spin_i_nd * res.n_i_loc
        M_ir_Y_add = -M_ir_spin_loc[2] * sθ - M_ir_spin_loc[3] * cθ
        M_ir_Z_add =  M_ir_spin_loc[2] * cθ - M_ir_spin_loc[3] * sθ

        M_ir_y += res.Z_rel_i * F_ir_X_add - res.X_rel_i * F_ir_Z_add + M_ir_Y_add
        M_ir_z += res.X_rel_i * F_ir_Y_add - res.Y_rel_i * F_ir_X_add + M_ir_Z_add

        # α2 states eliminated
    end
    # ── IR velocity derivatives (Newton-Euler) ──
    du.ir.vel.x  = F_ir_x / m_ir
    du.ir.vel.y  = F_ir_y / m_ir
    du.ir.vel.z  = F_ir_z / m_ir
    I_p_ir = m_ir * (D_i_star / 2)^2
    I_d_ir = I_p_ir / 2.0
    du.ir.vel.γy = (M_ir_y - I_p_ir * ctx.ω_ir * ctx.γ̇z_ir) / I_d_ir
    du.ir.vel.γz = (M_ir_z + I_p_ir * ctx.ω_ir * ctx.γ̇y_ir) / I_d_ir

    # ── Cage velocity derivatives ──
    du.cage.vel.x = F_cg_x / cage_mass
    du.cage.vel.y = F_cg_y / cage_mass
    du.cage.vel.z = F_cg_z / cage_mass
    du.cage.vel.θ = F_cg_θ / cage_Ixx

    # ── 热网络导数 ──
    T_amb = p.thermal.T_amb

    Q_ir_ball = p.thermal.G_ir_ball * (ctx.T_ir_node - ctx.T_ball_node)
    Q_or_ball = p.thermal.G_or_ball * (ctx.T_or_node - ctx.T_ball_node)
    Q_ball_oil = p.thermal.G_ball_oil * (ctx.T_ball_node - ctx.T_oil_node)
    Q_or_amb = p.thermal.G_or_amb * (ctx.T_or_node - T_amb)
    Q_oil_amb = p.thermal.G_oil_amb * (ctx.T_oil_node - T_amb)
    Q_oil_flow = p.thermal.oil_flow_mdot_cp * (ctx.T_oil_node - p.thermal.T_oil_inlet)

    dT_ir_dt = (H_ir_total - Q_ir_ball) / p.thermal.mcp_i
    dT_or_dt = (H_or_total - Q_or_ball - Q_or_amb) / p.thermal.mcp_o
    dT_ball_dt = (H_ball_total + Q_ir_ball + Q_or_ball - Q_ball_oil) / p.thermal.mcp_ball
    dT_oil_dt = (H_oil_total + Q_ball_oil - Q_oil_amb - Q_oil_flow) / p.thermal.mcp_oil

    du.thermal.T_i = dT_ir_dt
    du.thermal.T_o = dT_or_dt
    du.thermal.T_b = dT_ball_dt
    du.thermal.T_oil = dT_oil_dt

    # ── 热量积分器 ──
    du.heat.ir    = H_ir_total
    du.heat.or_   = H_or_total
    du.heat.b     = H_ball_total
    du.heat.oil   = H_oil_total
    du.heat.Q_amb = Q_or_amb + Q_oil_amb
    du.heat.Cp    = p.thermal.mcp_i * dT_ir_dt + p.thermal.mcp_o * dT_or_dt + p.thermal.mcp_ball * dT_ball_dt + p.thermal.mcp_oil * dT_oil_dt

    return nothing
end
