text = read("src/Dynamics/kernel.jl", String)
text = replace(text, 
    "        off_a2 = ball_alpha2_offset(j, Z)\n        du[off_a2]   = (res.alpha2_kin_i - u[off_a2]) / res.tau_a2_i\n        du[off_a2+1] = (res.alpha2_kin_o - u[off_a2+1]) / res.tau_a2_o" => "        # α2 states eliminated"
)
write("src/Dynamics/kernel.jl", text)
println("Patched kernel.jl")

text2 = read("src/Physics/kinematics.jl", String)
text2 = replace(text2, 
    "    off_a2_j = ball_alpha2_offset(j, Z)\n    a2_i_state = smooth_signed_saturate(u[off_a2_j], T_u(0.5))\n    a2_o_state = smooth_signed_saturate(u[off_a2_j + 1], T_u(0.5))" => "    a2_i_state = zero(T_u)\n    a2_o_state = zero(T_u)"
)
write("src/Physics/kinematics.jl", text2)
println("Patched kinematics.jl")

text3 = read("src/Dynamics/state.jl", String)
text3 = replace(text3, 
    "const N_ALPHA2_PER_BALL = 2" => "const N_ALPHA2_PER_BALL = 0"
)
text3 = replace(text3, 
    "         vel=(x=0.0, r=0.0, θ=0.0, ωx=0.0, ωy=0.0, ωz=0.0),\n         alpha2=(i=0.0, o=0.0))" => "         vel=(x=0.0, r=0.0, θ=0.0, ωx=0.0, ωy=0.0, ωz=0.0))"
)
text3 = replace(text3,
    "@inline ball_alpha2_view(u, j, Z) = u.ball[j].alpha2" => "# @inline ball_alpha2_view(u, j, Z) = u.ball[j].alpha2"
)
write("src/Dynamics/state.jl", text3)
println("Patched state.jl")
