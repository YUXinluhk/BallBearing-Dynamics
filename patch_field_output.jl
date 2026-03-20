lines = readlines("src/Dynamics/field_output.jl")

start_idx = findfirst(l -> occursin("── 内外圈状态", l), lines)
end_idx = findfirst(l -> occursin("── 写出场输出", l), lines) - 1

new_chunk = """
    ctx = build_kinematic_context(u, p, t)

    @inbounds for j in 1:Z
        base = (j - 1) * N_FIELD_PER_BALL
        res = compute_ball_tehd_kinematics(ctx, j)

        bp = ball_pos_view(u, j, Z)
        bv = ball_vel_view(u, j, Z)
        x_b, r_b, θ_b = bp[1], bp[2], bp[3]
        ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]
        q = ball_quat(u, j, Z)

        α_i = atan(res.dx_i, res.dr_i)
        α_o = atan(res.dx_o, res.dr_o)

        h_film_i = film_thickness_hd(res.u_entrain_i, res.Q_i_dim, ctx.E_prime_dim, res.R_ball_k * ctx.L_scale,
            (res.b_i_dim > 0 ? res.a_i_dim / res.b_i_dim : 1.0), ctx.μ_oil_dyn, p[P_ALPHA_PV], p[P_KTH], ctx.T_oil_node)
        h_film_o = film_thickness_hd(res.u_entrain_o, res.Q_o_dim, ctx.E_prime_dim, res.R_ball_k * ctx.L_scale,
            (res.b_o_dim > 0 ? res.a_o_dim / res.b_o_dim : 1.0), ctx.μ_oil_dyn, p[P_ALPHA_PV], p[P_KTH], ctx.T_oil_node)

        κ_i = res.F_trac_i_dim / (res.Q_i_dim + 1e-24)
        κ_o = res.F_trac_o_dim / (res.Q_o_dim + 1e-24)

        pocket_θ = ctx.θ_cg + (j - 1) * 2π / Z
        Δθ = θ_b - pocket_θ
        Δθ -= 2π * round(Δθ / (2π))
        gap = r_b * Δθ
        pen = hypot(gap, 1e-12) - p[P_POCKET_CLR]
        F_pk = pen > 0 ? p[P_K_POCKET] * pen * ctx.Q_scale : 0.0

        H_drag = (res.P_drag_nd + res.P_ehl_roll_nd) * ctx.Q_scale * ctx.V_scale
        H_churn = res.P_drag_spin_nd * ctx.Q_scale * ctx.V_scale
"""

open("src/Dynamics/field_output.jl", "w") do io
    for i in 1:(start_idx-1)
        println(io, lines[i])
    end
    print(io, new_chunk)
    for i in (end_idx+1):length(lines)
        line = lines[i]
        line = replace(line, "δ_i_sm * L_scale" => "res.δ_i_sm * ctx.L_scale")
        line = replace(line, "δ_o_sm * L_scale" => "res.δ_o_sm * ctx.L_scale")
        line = replace(line, "a_i_dim" => "res.a_i_dim")
        line = replace(line, "b_i_dim" => "res.b_i_dim")
        line = replace(line, "a_o_dim" => "res.a_o_dim")
        line = replace(line, "b_o_dim" => "res.b_o_dim")
        line = replace(line, "u_roll_i_dim" => "res.u_roll_i_dim")
        line = replace(line, "u_roll_o_dim" => "res.u_roll_o_dim")
        
        # Fixed exact replacements for the lines in field_output.jl
        line = replace(line, "u_mean_i * V_scale" => "0.5 * (res.v_ir_θ_dim + res.v_ball_θ * ctx.V_scale)")
        line = replace(line, "u_mean_o * V_scale" => "0.5 * (res.v_or_θ_dim + res.v_ball_θ * ctx.V_scale)")
        
        line = replace(line, "F_trac_i_dim" => "res.F_trac_i_dim")
        line = replace(line, "F_trac_o_dim" => "res.F_trac_o_dim")
        line = replace(line, "M_sp_i_dim" => "res.M_sp_i_dim")
        line = replace(line, "M_sp_o_dim" => "res.M_sp_o_dim")
        line = replace(line, "F_d_dim" => "res.F_drag_dim")
        line = replace(line, "M_ch_dim" => "res.M_ch_dim")
        
        line = replace(line, "ω_orb_x * W_scale" => "res.ω_ball_x_dim")
        line = replace(line, "ω_orb_r * W_scale" => "res.ω_ball_y_dim")
        line = replace(line, "ω_orb_θ * W_scale" => "res.ω_ball_z_dim")
        
        line = replace(line, "r_b * L_scale" => "r_b * ctx.L_scale")
        line = replace(line, "θ̇_b * W_scale" => "θ̇_b * ctx.W_scale")
        
        line = replace(line, "H_slide_i" => "res.H_slide_i_dim")
        line = replace(line, "H_slide_o" => "res.H_slide_o_dim")
        line = replace(line, "H_spin_i" => "res.H_spin_i_dim")
        line = replace(line, "H_spin_o" => "res.H_spin_o_dim")
        line = replace(line, "Q_i_dim" => "res.Q_i_dim")
        line = replace(line, "Q_o_dim" => "res.Q_o_dim")
        line = replace(line, "ω_spin_i_dim" => "res.ω_spin_i_dim")
        line = replace(line, "ω_spin_o_dim" => "res.ω_spin_o_dim")
        
        println(io, line)
    end
end
