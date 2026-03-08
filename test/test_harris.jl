# Test: 1D kernel model WITH Harris rolling speed
# = what the kernel SHOULD compute if ω_body had the Harris values

d_ball = 12.7e-3;
d_m = 70.0e-3;
α = deg2rad(40.0);
γ = d_ball / d_m;
n_ir = 10000.0;
ω_ir = n_ir * 2π / 60;

# Harris 
sα = sin(α);
cα = cos(α);
β = atan(sα, cα + γ)
sβ = sin(β);
cβ = cos(β);
tβ = sβ / cβ;
denom_R = ((cα + tβ * sα) / (1.0 + γ * cα) +
           (cα + tβ * sα) / (1.0 - γ * cα)) * γ * cβ
ω_R = -ω_ir / denom_R

# Nondim parameters (using L_scale = d_m for nondim)
L_scale = d_m
V_scale = L_scale * abs(ω_ir)  # V = L × W ≈ L × ω_ir
D_star = d_ball / L_scale
d_m_star = 1.0  # d_m / L_scale = 1
D_i_star = (d_m - d_ball * cα) / L_scale

# Nondim ω  
W_scale = abs(ω_ir)
ω_ir_nd = 1.0  # ω_ir / W_scale = 1

# Cage speed
ω_m = 0.5 * ω_ir * (1.0 - γ * cα)
ω_m_nd = ω_m / W_scale

# Ball center
r_b_nd = d_m / (2 * L_scale)  # = 0.5

# KERNEL 1D MODEL
R_ball_k = D_star / 2
v_ir_θ = (D_i_star / 2) * ω_ir_nd       # inner race tang vel (nondim)
v_ball_θ = r_b_nd * ω_m_nd              # ball orbital velocity (nondim)

# Harris ω in body frame (nondim)
wa0_nd = ω_R / W_scale * cβ    # ωx nondim
wa2_nd = ω_R / W_scale * sβ    # ωz nondim

# Rolling component at inner contact: ωx·cos(αi) - ωz·sin(αi) = ω_R/W × cos(αi - β)
ω_roll_i = wa0_nd * cα - wa2_nd * sα

# Ball surface velocity (1D kernel model)  
v_ball_surface_i = v_ball_θ - R_ball_k * ω_roll_i

# Inner race velocity at contact: v_race = (D_i/2) × ω_ir  
# But this is the nondim velocity (v_ir_θ already computed above)

u_slide_i_nd = abs(v_ir_θ - v_ball_surface_i)
u_slide_i_dim = u_slide_i_nd * V_scale

println("=== 1D Kernel Model with Harris ω ===")
println("v_ir_θ (nondim) = $(round(v_ir_θ, digits=6))")
println("v_ball_θ (nondim) = $(round(v_ball_θ, digits=6))")
println("R_ball_k = $(round(R_ball_k, digits=6))")
println("wa0_nd = $(round(wa0_nd, digits=4)), wa2_nd = $(round(wa2_nd, digits=4))")
println("ω_roll_i (nondim) = $(round(ω_roll_i, digits=4))")
println("v_ball_surface_i (nondim) = $(round(v_ball_surface_i, digits=6))")
println("u_slide_i (nondim) = $(round(u_slide_i_nd, digits=6))")
println("u_slide_i (dim) = $(round(u_slide_i_dim, digits=4)) m/s")

# Now also test with original ODE ω (wa0 ≈ 0, wa2 ≈ orbital)
wa0_ode = 0.0
wa2_ode = ω_m_nd
ω_roll_i_ode = wa0_ode * cα - wa2_ode * sα
v_ball_surface_i_ode = v_ball_θ - R_ball_k * ω_roll_i_ode
u_slide_i_ode = abs(v_ir_θ - v_ball_surface_i_ode) * V_scale
println("\n--- ODE original (wa0≈0, wa2≈ω_m) ---")
println("ω_roll_i_ode = $(round(ω_roll_i_ode, digits=4))")
println("v_surf_i_ode = $(round(v_ball_surface_i_ode, digits=6))")
println("u_slide_i_ode = $(round(u_slide_i_ode, digits=4)) m/s")

# What ω_roll_i would give ZERO slide?
# v_ir_θ = v_ball_θ - R × ω_roll_pure  
# ω_roll_pure = (v_ball_θ - v_ir_θ) / R
ω_roll_pure = (v_ball_θ - v_ir_θ) / R_ball_k
println("\n--- Pure rolling target ---")
println("ω_roll for zero slide = $(round(ω_roll_pure, digits=4))")
println("ω_R for zero slide = $(round(ω_roll_pure * W_scale / cos(α - β), digits=1)) rad/s")
println("Harris ω_R = $(round(ω_R, digits=1)) rad/s")
