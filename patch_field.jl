text = read("src/Dynamics/field_output.jl", String)

old_setup = \"\"\"
    # ── 内外圈状态 (完全镜像 kernel.jl) ──
    T_u = Float64
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
\"\"\"
new_setup = \"\"\"
    # ── 内外圈状态 ──
    T_u = Float64
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
        q = ball_quat(u, j, Z)
\"\"\"
new_ball = \"\"\"
        bp = u.ball[j].pos
        bv = u.ball[j].vel
        x_b, r_b, θ_b = bp.x, bp.r, bp.θ
        ẋ_b, ṙ_b, θ̇_b = bv.x, bv.r, bv.θ
        q = ball_quat(u, j, Z)
\"\"\"
text = replace(text, old_ball => new_ball)

open("src/Dynamics/field_output.jl", "w") do io
    write(io, text)
end
println("Patched field_output.jl!")
