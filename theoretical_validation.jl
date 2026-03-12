# theoretical_validation.jl
using Pkg
Pkg.activate(@__DIR__)
using TOML
using Printf

# Include core parser
include(joinpath(@__DIR__, "src", "ADORE.jl"))
using .ADORE

if length(ARGS) < 1
    println("Usage: julia theoretical_validation.jl <case.toml>")
    exit(1)
end

case_file = ARGS[1]
config_file = joinpath(@__DIR__, "inputs", case_file)

if !isfile(config_file)
    println("Error: File $config_file not found.")
    exit(1)
end

println("="^60)
println("  ADORE V2 -> Theoretical Validation vs. Classical Mechanics")
println("="^60)
println("Analyzing Case: $case_file\n")

# Load configuration using the new TOML parser
geom, mat, lub, trac, cage, config = load_simulation_config(config_file)

# --- 1. Geometry & Kinematics (Harris Theory) ---
d_m = geom.d_m        # Pitch diameter
D_w = geom.d          # Ball diameter
f_i = geom.f_i        # Inner race curvature
f_o = geom.f_o        # Outer race curvature
alpha = geom.alpha_0  # Nominal contact angle

n_rpm = config.inner_race_speed * 30 / π
omega_i = n_rpm * π / 30

gamma = D_w * cos(alpha) / d_m

# Theoretical Cage Speed (ω_c)
omega_c_theory = 0.5 * omega_i * (1 - gamma)
n_c_theory = omega_c_theory * 30 / π

# Theoretical Ball Orbital Speed (ω_m) = Cage Speed
# Theoretical Ball Spin Speed (ω_B)
omega_B_theory = 0.5 * omega_i * (d_m / D_w) * (1 - gamma^2)
n_B_theory = omega_B_theory * 30 / π

@printf("--- KINEMATIC VALIDATION (Harris) ---\n")
@printf("Input Inner Race Speed : %.2f RPM\n", n_rpm)
@printf("Theoretical Cage Speed : %.2f RPM (%.4f rad/s)\n", n_c_theory, omega_c_theory)
@printf("Theoretical Ball Spin  : %.2f RPM (%.4f rad/s)\n", n_B_theory, omega_B_theory)
@printf("Pitch over Ball Ratio (d_m/D_w) : %.2f\n", d_m / D_w)
println()

# --- 2. Load & Hertzian Contact (Simple Approximation) ---
F_a = config.F_axial
F_r = config.F_radial
Z = geom.n_balls

@printf("--- LOAD DISTRIBUTION VALIDATION ---\n")
@printf("Applied Axial Load     : %.2f N\n", F_a)
@printf("Applied Radial Load    : %.2f N\n", F_r)

if F_r == 0 && F_a > 0
    # Pure axial load - load is distributed equally
    Q_avg = F_a / (Z * sin(alpha))
    @printf("Theoretical Ball Load (Pure Axial) : %.2f N\n", Q_avg)
elseif F_a == 0 && F_r > 0
    # Pure radial load
    # Stribeck's constant approximation: Q_max = 4.37 * F_r / (Z * cos(alpha))
    # Or more commonly: 5 * F_r / Z
    Q_max = 5.0 * F_r / (Z * cos(alpha))
    @printf("Max Theoretical Ball Load (Stribeck) : ~%.2f N\n", Q_max)
else
    @printf("Combined Load - Requires numerical iteration for exact max load.\n")
end
println()

# --- 3. Friction Torque (Palmgren Empirical) ---
# For angular contact ball bearings
# Kinematic viscosity ν [cSt = mm^2/s] = Dynamic viscosity [Pa·s = kg/(m·s)] / Density [kg/m^3] * 1e6
nu = (lub.mu_0 / lub.rho_lub) * 1e6 # cSt (mm^2/s)
d_m_mm = d_m * 1000

# 1. Viscous Friction Torque (M_v)
f_0 = 2.0 # factor for angular contact
if nu * n_rpm > 2000
    M_v = 1e-7 * f_0 * (nu * n_rpm)^(2 / 3) * d_m_mm^3 # N-mm
else
    M_v = 160 * 1e-7 * f_0 * d_m_mm^3
end
M_v /= 1000 # Convert to N-m

# 2. Load Friction Torque (M_l)
f_1 = 0.001 # approx for angular contact
P_1 = max(F_r, F_a) # Simplistic equivalent load
M_l = f_1 * P_1 * d_m_mm
M_l /= 1000 # Convert to N-m

M_tot = M_v + M_l

# Power Loss
Power_loss = M_tot * omega_i

@printf("--- LOW-SPEED FRICTION & POWER LOSS (Palmgren) ---\n")
@printf("Kinematic Viscosity (ν) : %.2f cSt\n", nu)
@printf("Viscous Torque (M_v)    : %.4f N-m\n", M_v)
@printf("Load Torque (M_l)       : %.4f N-m\n", M_l)
@printf("Total Frictional Torque : %.4f N-m\n", M_tot)
@printf("Total Heat Generation   : %.2f W\n", Power_loss)
println()

println("=> To compare these with the numerical ADORE integration, check the final console output of 'run_single.jl' or the generated plots.")
