text = read("src/Dynamics/kernel.jl", String)

# the view setup
old_setup = \"\"\"
    ir_p = ir_pos_view(u, Z)
    ir_v = ir_vel_view(u, Z)
    x_ir, y_ir, z_ir, γy_ir, γz_ir = ir_p[1], ir_p[2], ir_p[3], ir_p[4], ir_p[5]
    ẋ_ir, ẏ_ir, ż_ir, γ̇y_ir, γ̇z_ir = ir_v[1], ir_v[2], ir_v[3], ir_v[4], ir_v[5]

    cg_v = cage_vel_view(u, Z)
    θ̇_cg = cg_v[4]
    cg_p = cage_pos_view(u, Z)
    θ_cg = cg_p[4]

    # Thermal node temperatures
    tv = thermal_view(u, Z)
    T_ir_t, T_or_t = tv[1], tv[2]
    T_ball_t = tv[3]
    T_oil_t = tv[4]
\"\"\"
new_setup = \"\"\"
    ir_p = u.ir.pos
    ir_v = u.ir.vel
    x_ir, y_ir, z_ir, γy_ir, γz_ir = ir_p.x, ir_p.y, ir_p.z, ir_p.γy, ir_p.γz
    ẋ_ir, ẏ_ir, ż_ir, γ̇y_ir, γ̇z_ir = ir_v.x, ir_v.y, ir_v.z, ir_v.γy, ir_v.γz

    cg_v = u.cage.vel
    θ̇_cg = cg_v.θ
    cg_p = u.cage.pos
    θ_cg = cg_p.θ

    # Thermal node temperatures
    tv = u.thermal
    T_ir_t, T_or_t = tv.Ti, tv.To
    T_ball_t = tv.Tb
    T_oil_t = tv.Toil
\"\"\"
text = replace(text, old_setup => new_setup)

old_ball = \"\"\"
        bp = ball_pos_view(u, j, Z)
        bv = ball_vel_view(u, j, Z)
        x_b, r_b, θ_b = bp[1], bp[2], bp[3]
        ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]
        ω_orb = ball_omega(u, j, Z)
        q = ball_quat(u, j, Z)
\"\"\"
new_ball = \"\"\"
        bp = u.ball[j].pos
        bv = u.ball[j].vel
        x_b, r_b, θ_b = bp.x, bp.r, bp.θ
        ẋ_b, ṙ_b, θ̇_b = bv.x, bv.r, bv.θ
        ω_orb = ball_omega(u, j, Z)
        q = ball_quat(u, j, Z)
\"\"\"
text = replace(text, old_ball => new_ball)

old_a2 = \"\"\"
        off_a2_j = ball_alpha2_offset(j, Z)
        a2_i_state = smooth_signed_saturate(u[off_a2_j], 0.5)
        a2_o_state = smooth_signed_saturate(u[off_a2_j + 1], 0.5)
\"\"\"
new_a2 = \"\"\"
        a2_i_state = smooth_signed_saturate(u.ball[j].alpha2.i, 0.5)
        a2_o_state = smooth_signed_saturate(u.ball[j].alpha2.o, 0.5)
\"\"\"
text = replace(text, old_a2 => new_a2)

# Write back du for balls
old_du_ball = \"\"\"
        # ---- Write derivatives ----
        du_bp = ball_pos_view(du, j, Z)
        du_bp[1] = ẋ_b * (ctx.T_scale / ctx.L_scale)
        du_bp[2] = ṙ_b * (ctx.T_scale / ctx.L_scale)
        du_bp[3] = θ̇_b * ctx.T_scale

        q_dot_scalar = 0.5 * (-q.v1 * ω_b_b[1] - q.v2 * ω_b_b[2] - q.v3 * ω_b_b[3])
        q_dot_v1 = 0.5 * (q.s * ω_b_b[1] + q.v2 * ω_b_b[3] - q.v3 * ω_b_b[2])
        q_dot_v2 = 0.5 * (q.s * ω_b_b[2] + q.v3 * ω_b_b[1] - q.v1 * ω_b_b[3])
        q_dot_v3 = 0.5 * (q.s * ω_b_b[3] + q.v1 * ω_b_b[2] - q.v2 * ω_b_b[1])

        du_bp[4] = q_dot_scalar * ctx.T_scale
        du_bp[5] = q_dot_v1 * ctx.T_scale
        du_bp[6] = q_dot_v2 * ctx.T_scale
        du_bp[7] = q_dot_v3 * ctx.T_scale

        du_bv = ball_vel_view(du, j, Z)
        du_bv[1] = ax_b_dim * (ctx.T_scale^2 / ctx.L_scale)
        du_bv[2] = ar_b_dim * (ctx.T_scale^2 / ctx.L_scale)
        du_bv[3] = aθ_b_dim * ctx.T_scale^2

        du_bv[4] = res.ω_dot_x_dim * ctx.T_scale^2
        du_bv[5] = res.ω_dot_y_dim * ctx.T_scale^2
        du_bv[6] = res.ω_dot_z_dim * ctx.T_scale^2

        du[off_a2_j] = (res.α2_kin_i - u[off_a2_j]) / (tau_a2_star_i + 1e-12)
        du[off_a2_j + 1] = (res.α2_kin_o - u[off_a2_j + 1]) / (tau_a2_star_o + 1e-12)
\"\"\"
new_du_ball = \"\"\"
        # ---- Write derivatives ----
        du_bp = du.ball[j].pos
        du_bp.x = ẋ_b * (ctx.T_scale / ctx.L_scale)
        du_bp.r = ṙ_b * (ctx.T_scale / ctx.L_scale)
        du_bp.θ = θ̇_b * ctx.T_scale

        q_dot_scalar = 0.5 * (-q.v1 * ω_b_b[1] - q.v2 * ω_b_b[2] - q.v3 * ω_b_b[3])
        q_dot_v1 = 0.5 * (q.s * ω_b_b[1] + q.v2 * ω_b_b[3] - q.v3 * ω_b_b[2])
        q_dot_v2 = 0.5 * (q.s * ω_b_b[2] + q.v3 * ω_b_b[1] - q.v1 * ω_b_b[3])
        q_dot_v3 = 0.5 * (q.s * ω_b_b[3] + q.v1 * ω_b_b[2] - q.v2 * ω_b_b[1])

        du_bp.q0 = q_dot_scalar * ctx.T_scale
        du_bp.q1 = q_dot_v1 * ctx.T_scale
        du_bp.q2 = q_dot_v2 * ctx.T_scale
        du_bp.q3 = q_dot_v3 * ctx.T_scale

        du_bv = du.ball[j].vel
        du_bv.x = ax_b_dim * (ctx.T_scale^2 / ctx.L_scale)
        du_bv.r = ar_b_dim * (ctx.T_scale^2 / ctx.L_scale)
        du_bv.θ = aθ_b_dim * ctx.T_scale^2

        du_bv.ωx = res.ω_dot_x_dim * ctx.T_scale^2
        du_bv.ωy = res.ω_dot_y_dim * ctx.T_scale^2
        du_bv.ωz = res.ω_dot_z_dim * ctx.T_scale^2

        du.ball[j].alpha2.i = (res.α2_kin_i - u.ball[j].alpha2.i) / (tau_a2_star_i + 1e-12)
        du.ball[j].alpha2.o = (res.α2_kin_o - u.ball[j].alpha2.o) / (tau_a2_star_o + 1e-12)
\"\"\"
text = replace(text, old_du_ball => new_du_ball)

# Write back du for IR, cage, thermal
old_du_ir = \"\"\"
    du_ir_p = ir_pos_view(du, Z)
    du_ir_v = ir_vel_view(du, Z)
    du_ir_p[1] = ẋ_ir_dim * ctx.T_scale / ctx.L_scale
    du_ir_p[2] = ẏ_ir_dim * ctx.T_scale / ctx.L_scale
    du_ir_p[3] = ż_ir_dim * ctx.T_scale / ctx.L_scale
    du_ir_p[4] = γ̇y_ir_dim * ctx.T_scale
    du_ir_p[5] = γ̇z_ir_dim * ctx.T_scale

    du_ir_v[1] = x_ddot_ir * ctx.T_scale^2 / ctx.L_scale
    du_ir_v[2] = y_ddot_ir * ctx.T_scale^2 / ctx.L_scale
    du_ir_v[3] = z_ddot_ir * ctx.T_scale^2 / ctx.L_scale
    du_ir_v[4] = γy_ddot_ir * ctx.T_scale^2
    du_ir_v[5] = γz_ddot_ir * ctx.T_scale^2

    du_cg_p = cage_pos_view(du, Z)
    du_cg_v = cage_vel_view(du, Z)
    du_cg_p[1] = ẋ_cg_dim * ctx.T_scale / ctx.L_scale
    du_cg_p[2] = ẏ_cg_dim * ctx.T_scale / ctx.L_scale
    du_cg_p[3] = ż_cg_dim * ctx.T_scale / ctx.L_scale
    du_cg_p[4] = θ̇_cg_dim * ctx.T_scale

    du_cg_v[1] = x_ddot_cg * ctx.T_scale^2 / ctx.L_scale
    du_cg_v[2] = y_ddot_cg * ctx.T_scale^2 / ctx.L_scale
    du_cg_v[3] = z_ddot_cg * ctx.T_scale^2 / ctx.L_scale
    du_cg_v[4] = θ_ddot_cg * ctx.T_scale^2

    du_tv = thermal_view(du, Z)
    # thermal ODE returns [K/s], scale to nondim
    du_tv .= dT_dt .* ctx.T_scale

    # Heat accumulators
    base_ha = heat_accum_offset(Z)
    du[base_ha] = F_drag_total_dim * (ctx.T_scale / ctx.Q_scale)
    du[base_ha+1] = 0.0
    du[base_ha+2] = sum_H_spin_i_dim * (ctx.T_scale / (ctx.Q_scale * ctx.V_scale))
    du[base_ha+3] = H_oil_total * (ctx.T_scale / (ctx.Q_scale * ctx.V_scale))
    du[base_ha+4] = 0.0 # etc
    du[base_ha+5] = 0.0
\"\"\"

new_du_ir = \"\"\"
    du.ir.pos.x = ẋ_ir_dim * ctx.T_scale / ctx.L_scale
    du.ir.pos.y = ẏ_ir_dim * ctx.T_scale / ctx.L_scale
    du.ir.pos.z = ż_ir_dim * ctx.T_scale / ctx.L_scale
    du.ir.pos.γy = γ̇y_ir_dim * ctx.T_scale
    du.ir.pos.γz = γ̇z_ir_dim * ctx.T_scale

    du.ir.vel.x = x_ddot_ir * ctx.T_scale^2 / ctx.L_scale
    du.ir.vel.y = y_ddot_ir * ctx.T_scale^2 / ctx.L_scale
    du.ir.vel.z = z_ddot_ir * ctx.T_scale^2 / ctx.L_scale
    du.ir.vel.γy = γy_ddot_ir * ctx.T_scale^2
    du.ir.vel.γz = γz_ddot_ir * ctx.T_scale^2

    du.cage.pos.x = ẋ_cg_dim * ctx.T_scale / ctx.L_scale
    du.cage.pos.y = ẏ_cg_dim * ctx.T_scale / ctx.L_scale
    du.cage.pos.z = ż_cg_dim * ctx.T_scale / ctx.L_scale
    du.cage.pos.θ = θ̇_cg_dim * ctx.T_scale

    du.cage.vel.x = x_ddot_cg * ctx.T_scale^2 / ctx.L_scale
    du.cage.vel.y = y_ddot_cg * ctx.T_scale^2 / ctx.L_scale
    du.cage.vel.z = z_ddot_cg * ctx.T_scale^2 / ctx.L_scale
    du.cage.vel.θ = θ_ddot_cg * ctx.T_scale^2

    du.thermal .= dT_dt .* ctx.T_scale

    # Heat accumulators
    du.heat.ir = F_drag_total_dim * (ctx.T_scale / ctx.Q_scale)
    du.heat.or_ = 0.0
    du.heat.b = sum_H_spin_i_dim * (ctx.T_scale / (ctx.Q_scale * ctx.V_scale))
    du.heat.oil = H_oil_total * (ctx.T_scale / (ctx.Q_scale * ctx.V_scale))
    du.heat.Q_amb = 0.0 # etc
    du.heat.Cp_dT = 0.0
\"\"\"
text = replace(text, old_du_ir => new_du_ir)

open("src/Dynamics/kernel.jl", "w") do io
    write(io, text)
end
println("Patched kernel.jl!")
