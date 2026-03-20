text = read("src/Dynamics/state.jl", String)

# the views section
start_v = findfirst(l -> occursin("── Zero-copy views", l), split(text, '\n'))
end_v = findfirst(l -> occursin("── State initialization", l), split(text, '\n')) - 1

lines = split(text, '\n')
old_views = join(lines[start_v:end_v], '\n') * '\n'

new_views = """
# ── ComponentArrays Properties ───────────────────────────────────────

"View inner race position [x, y, z, γ_y, γ_z]"
@inline ir_pos_view(u, Z) = u.ir.pos

"View ball j position [x, r, θ, q₀, q₁, q₂, q₃]"
@inline ball_pos_view(u, j, Z) = u.ball[j].pos

"View ball j translational position [x, r, θ]"
@inline ball_trans_view(u, j, Z) = @view u.ball[j].pos[1:3]

"Extract ball j quaternion as Quaternion"
@inline function ball_quat(u, j, Z)
    bp = u.ball[j].pos
    Quaternion{typeof(bp[1])}(bp.q0, bp.q1, bp.q2, bp.q3)
end

"View cage position [x, y, z, θ_cage]"
@inline cage_pos_view(u, Z) = u.cage.pos

"View inner race velocity [ẋ, ẏ, ż, γ̇_y, γ̇_z]"
@inline ir_vel_view(u, Z) = u.ir.vel

"View ball j velocity [ẋ, ṙ, θ̇, ω_x, ω_y, ω_z]"
@inline ball_vel_view(u, j, Z) = u.ball[j].vel

"Ball j angular velocity [ω_x, ω_y, ω_z] as SVector"
@inline function ball_omega(u, j, Z)
    bv = u.ball[j].vel
    SVector{3,typeof(bv[1])}(bv.ωx, bv.ωy, bv.ωz)
end

"View cage velocity [ẋ, ẏ, ż, θ̇_cage]"
@inline cage_vel_view(u, Z) = u.cage.vel

"View ball j internal α₂ states [α₂_i, α₂_o]"
@inline ball_alpha2_view(u, j, Z) = u.ball[j].alpha2

"View thermal nodes [T_i, T_o, T_b, T_oil]"
@inline thermal_view(u, Z) = u.thermal

# ── Quaternion write-back ────────────────────────────────────────────

"Set ball j quaternion in state vector"
@inline function set_ball_quat!(u, j, Z, q::Quaternion)
    u.ball[j].pos.q0 = real(q)
    i1, i2, i3 = imag_part(q)
    u.ball[j].pos.q1 = i1
    u.ball[j].pos.q2 = i2
    u.ball[j].pos.q3 = i3
    return nothing
end

"""

text = replace(text, old_views => new_views)

start_i = findfirst(l -> occursin("init_state(geom::BearingGeometry", l), split(text, '\n'))
old_init = join(split(text, '\n')[start_i:end], '\n')

new_init = """
init_state(geom::BearingGeometry, qs::QuasiStaticResult,
    config::SimulationConfig, lub::LubricantParams)
    Z = geom.n_balls

    ball_states = map(1:Z) do j
        (pos=(x=0.0, r=0.0, θ=0.0, q0=1.0, q1=0.0, q2=0.0, q3=0.0),
         vel=(x=0.0, r=0.0, θ=0.0, ωx=0.0, ωy=0.0, ωz=0.0),
         alpha2=(i=0.0, o=0.0))
    end
    
    u0_namedtuple = (
        ir = (pos=(x=qs.delta_a, y=qs.delta_ry, z=qs.delta_rz, γy=0.0, γz=0.0), vel=(x=0.0, y=0.0, z=0.0, γy=0.0, γz=0.0)),
        cage = (pos=(x=0.0, y=0.0, z=0.0, θ=0.0), vel=(x=0.0, y=0.0, z=0.0, θ=mean(qs.omega_m))),
        ball = ball_states,
        thermal = (T_i=config.thermal.T_init, T_o=config.thermal.T_init, T_b=config.thermal.T_init, T_oil=config.thermal.T_init),
        heat = (ir=0.0, or_=0.0, b=0.0, oil=0.0, Q_amb=0.0, Cp=0.0)
    )
    
    u0 = ComponentArray(u0_namedtuple)

    D = geom.d
    f_i = geom.f_i
    f_o = geom.f_o
    α₀ = alpha_free(geom)
    B_i = (f_i - 0.5) * D
    B_o = (f_o - 0.5) * D

    x_ci = B_i * sin(α₀)
    r_ci = geom.d_m / 2 + B_i * cos(α₀)
    x_co = -B_o * sin(α₀)
    r_co = geom.d_m / 2 - B_o * cos(α₀)

    for j in 1:Z
        ψ = (j - 1) * ball_spacing(geom)
        sψ, cψ = sincos(ψ)

        α_o_j = deg2rad(qs.alpha_outer[j])
        L_o_j = B_o + qs.delta_outer[j]
        x_from_o = x_co + L_o_j * sin(α_o_j)
        r_from_o = r_co + L_o_j * cos(α_o_j)

        α_i_j = deg2rad(qs.alpha_inner[j])
        L_i_j = B_i + qs.delta_inner[j]
        r_ci_local = r_ci - qs.delta_ry * sψ + qs.delta_rz * cψ
        x_from_i = (x_ci + qs.delta_a) - L_i_j * sin(α_i_j)
        r_from_i = r_ci_local - L_i_j * cos(α_i_j)

        bp = u0.ball[j].pos
        bp.x = 0.5 * (x_from_o + x_from_i)
        bp.r = 0.5 * (r_from_o + r_from_i)
        bp.θ = (j - 1) * ball_spacing(geom)
        
        bv = u0.ball[j].vel
        bv.θ = qs.omega_m[j]
        beta_j = deg2rad(qs.beta[j])
        bv.ωx = qs.omega_R[j] * cos(beta_j)
        bv.ωy = qs.omega_R[j] * sin(beta_j)
        bv.ωz = 0.0
    end

    return u0
end
"""
text = replace(text, old_init => new_init)

open("src/Dynamics/state.jl", "w") do io
    write(io, text)
end
println("Patched state.jl successfully via patch_state_v2.jl!")
