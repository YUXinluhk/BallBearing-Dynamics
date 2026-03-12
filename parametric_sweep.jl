# Multi-Speed Parametric Sweep for ADORE V2
# Generates cage speed, friction torque, and heat generation vs shaft speed
# For journal paper experimental validation section

using Pkg
Pkg.activate(".")

include("src/ADORE.jl")
using .ADORE
using Printf
using Plots
using TOML

# ============================================================
# Sweep Parameters
# ============================================================
speeds_rpm = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
    12000, 14000, 16000, 18000, 20000]
F_a = 1000.0  # Axial load [N]
F_r = 0.0     # Radial load [N]

# Results storage
n_speeds = length(speeds_rpm)
cage_speeds_rpm = zeros(n_speeds)
harris_cage_rpm = zeros(n_speeds)
total_heat_W = zeros(n_speeds)
friction_torque_Nm = zeros(n_speeds)

# Load base config
base_toml = TOML.parsefile("Case3.toml")

println("="^60)
println("  ADORE V2 — Multi-Speed Parametric Sweep")
println("  F_a = $F_a N, F_r = $F_r N")
println("="^60)

for (idx, n_rpm) in enumerate(speeds_rpm)
    println("\n--- Speed: $n_rpm RPM ($(idx)/$(n_speeds)) ---")

    if n_rpm == 0
        cage_speeds_rpm[idx] = 0.0
        harris_cage_rpm[idx] = 0.0
        total_heat_W[idx] = 0.0
        friction_torque_Nm[idx] = 0.0
        println("  Skipping zero-speed case (static)")
        continue
    end

    # Create modified config
    cfg = deepcopy(base_toml)
    cfg["operating"]["shaft_speed_rpm"] = n_rpm
    cfg["operating"]["F_axial"] = F_a
    cfg["operating"]["F_radial"] = F_r

    # Write temporary TOML
    temp_toml = joinpath("cases", "sweep_n$(n_rpm).toml")
    open(temp_toml, "w") do f
        TOML.print(f, cfg)
    end

    # Load config and run simulation
    try
        config = load_simulation_config(temp_toml)

        # Quasi-static solution
        qs = solve_quasi_static(config)

        # Harris theoretical cage speed
        d = config.geometry.d
        d_m = config.geometry.d_m
        alpha0 = config.geometry.alpha_0
        harris_cage_rpm[idx] = 0.5 * n_rpm * (1.0 - d * cosd(alpha0) / d_m)

        # Dynamic solution
        sol, p_vec, scales = solve_dynamics(config, qs)

        if sol.retcode == :Success || sol.retcode == :Terminated
            # Extract cage speed from final state
            Z = config.geometry.Z
            cv = cage_vel_view(sol.u[end], Z)
            omega_cage_nd = cv[4]  # nondim cage angular velocity
            W_scale = scales.W
            cage_speeds_rpm[idx] = omega_cage_nd * W_scale * 60.0 / (2π)

            # Extract heat generation from field output
            fo = compute_field_output(sol, p_vec, scales, config)
            # Total heat at last time step
            total_heat_W[idx] = sum(fo.H_total[end, :])

            # Friction torque from total heat / shaft speed
            omega_shaft = n_rpm * 2π / 60.0
            friction_torque_Nm[idx] = total_heat_W[idx] / omega_shaft

            @printf("  Cage: %.1f RPM (Harris: %.1f), Heat: %.2f W, Torque: %.4f Nm\n",
                cage_speeds_rpm[idx], harris_cage_rpm[idx], total_heat_W[idx], friction_torque_Nm[idx])
        else
            println("  ⚠ Solver did not converge: $(sol.retcode)")
        end
    catch e
        println("  ✗ Error at $n_rpm RPM: $e")
    end
end

# ============================================================
# Palmgren Empirical Model for Comparison
# ============================================================
d = base_toml["geometry"]["d"]
d_m = base_toml["geometry"]["d_m"]
alpha0 = base_toml["geometry"]["alpha_0"]
nu_oil = 12.0  # kinematic viscosity at operating temp [mm²/s] (typical)

palmgren_torque = zeros(n_speeds)
for (idx, n_rpm) in enumerate(speeds_rpm)
    if n_rpm == 0
        continue
    end
    # Palmgren friction torque model
    f0 = 2.0  # for ball bearings with oil lubrication
    f1 = 0.0004  # for angular contact ball bearings
    P1 = max(F_a * 3.0 * cosd(alpha0), F_a)  # equivalent load

    # M_v: viscous friction torque
    M_v = if nu_oil * n_rpm >= 2000
        1e-7 * f0 * (nu_oil * n_rpm)^(2 / 3) * (d_m * 1000)^3 * 1e-9  # convert mm to m
    else
        160e-7 * f0 * (d_m * 1000)^3 * 1e-9
    end

    # M_l: load-dependent friction torque  
    M_l = f1 * P1 * d_m

    palmgren_torque[idx] = M_v + M_l
end

# ============================================================
# SKF Model Approximation
# ============================================================
skf_torque = zeros(n_speeds)
for (idx, n_rpm) in enumerate(speeds_rpm)
    if n_rpm == 0
        continue
    end
    omega = n_rpm * 2π / 60.0
    # SKF new model: M = M_rr + M_sl + M_seal + M_drag
    # Simplified version for unsealed bearing

    # Rolling friction
    phi_ish = 1.0 / (1.0 + 1.84e-9 * (n_rpm * d_m * 1000)^1.28 * nu_oil^0.64)
    phi_rs = 1.0  # no seal
    G_rr = 3.73 * (F_a * d_m * 1000)^0.57  # for angular contact
    M_rr = phi_ish * phi_rs * G_rr * 1e-6

    # Sliding friction
    mu_sl = 0.05  # typical for grease/oil
    G_sl = 2.84 * (F_a * d_m * 1000)^0.54
    M_sl = G_sl * mu_sl * 1e-6

    # Drag loss
    V_m = 0.6e-6 * pi * d_m^3  # oil volume approximation
    M_drag = 1.093e-7 * n_rpm^2 * (d_m * 1000)^3 * (nu_oil / (nu_oil + 1e6))^0.5 * 1e-9

    skf_torque[idx] = M_rr + M_sl + M_drag
end

# ============================================================
# Generate Plots
# ============================================================
out_dir = "paper/figures"
mkpath(out_dir)

# Plot 1: Cage Speed vs Shaft Speed
fig1 = plot(speeds_rpm, harris_cage_rpm,
    label="Harris Theory", linestyle=:dash, linewidth=2, color=:black,
    xlabel="Shaft Speed (RPM)", ylabel="Cage Speed (RPM)",
    title="Cage Orbital Speed: ADORE vs Harris Theory",
    legend=:topleft, grid=true, gridalpha=0.3, size=(800, 500))
scatter!(fig1, speeds_rpm[speeds_rpm.>0], cage_speeds_rpm[speeds_rpm.>0],
    label="ADORE V2", markersize=6, color=:royalblue)
savefig(fig1, joinpath(out_dir, "cage_speed_vs_shaft.png"))
println("\n[OK] cage_speed_vs_shaft.png")

# Plot 2: Heat Generation vs Shaft Speed
fig2 = plot(speeds_rpm[speeds_rpm.>0], total_heat_W[speeds_rpm.>0],
    label="ADORE V2 (TEHD)", linewidth=2, marker=:circle, color=:red,
    xlabel="Shaft Speed (RPM)", ylabel="Total Heat Generation (W)",
    title="Bearing Heat Generation vs Shaft Speed",
    legend=:topleft, grid=true, gridalpha=0.3, size=(800, 500))
savefig(fig2, joinpath(out_dir, "heat_vs_speed.png"))
println("[OK] heat_vs_speed.png")

# Plot 3: Friction Torque Comparison
fig3 = plot(speeds_rpm[speeds_rpm.>0], friction_torque_Nm[speeds_rpm.>0] .* 1000,
    label="ADORE V2 (TEHD)", linewidth=2, marker=:circle, color=:royalblue,
    xlabel="Shaft Speed (RPM)", ylabel="Friction Torque (mN·m)",
    title="Friction Torque: ADORE vs Palmgren vs SKF",
    legend=:topleft, grid=true, gridalpha=0.3, size=(800, 500))
plot!(fig3, speeds_rpm[speeds_rpm.>0], palmgren_torque[speeds_rpm.>0] .* 1000,
    label="Palmgren (1959)", linestyle=:dash, linewidth=2, color=:darkgreen)
plot!(fig3, speeds_rpm[speeds_rpm.>0], skf_torque[speeds_rpm.>0] .* 1000,
    label="SKF Model (approx.)", linestyle=:dashdot, linewidth=2, color=:darkorange)
savefig(fig3, joinpath(out_dir, "friction_torque_comparison.png"))
println("[OK] friction_torque_comparison.png")

# Plot 4: Cage Speed Error vs Harris
cage_err = abs.(cage_speeds_rpm .- harris_cage_rpm) ./ max.(harris_cage_rpm, 1.0) .* 100
fig4 = plot(speeds_rpm[speeds_rpm.>0], cage_err[speeds_rpm.>0],
    label="ADORE vs Harris Error", linewidth=2, marker=:diamond, color=:purple,
    xlabel="Shaft Speed (RPM)", ylabel="Relative Error (%)",
    title="Cage Speed Prediction Error vs Harris Theory",
    legend=:topright, grid=true, gridalpha=0.3, size=(800, 500),
    ylims=(0, 2))
hline!(fig4, [0.5], linestyle=:dash, color=:gray, label="0.5% threshold", alpha=0.6)
savefig(fig4, joinpath(out_dir, "cage_speed_error.png"))
println("[OK] cage_speed_error.png")

# ============================================================
# Print Summary Table
# ============================================================
println("\n" * "="^80)
println("  PARAMETRIC SWEEP RESULTS SUMMARY")
println("="^80)
@printf("  %-8s  %-10s  %-10s  %-8s  %-10s  %-10s  %-10s\n",
    "RPM", "CageADORE", "CageHarris", "Err(%)", "Heat(W)", "M_ADORE", "M_Palm")
println("-"^80)
for idx in 1:n_speeds
    if speeds_rpm[idx] > 0
        err = abs(cage_speeds_rpm[idx] - harris_cage_rpm[idx]) / harris_cage_rpm[idx] * 100
        @printf("  %-8d  %-10.1f  %-10.1f  %-8.3f  %-10.2f  %-10.4f  %-10.4f\n",
            speeds_rpm[idx], cage_speeds_rpm[idx], harris_cage_rpm[idx], err,
            total_heat_W[idx], friction_torque_Nm[idx], palmgren_torque[idx])
    end
end
println("="^80)

# Save data to CSV for paper
open(joinpath(out_dir, "sweep_data.csv"), "w") do f
    println(f, "rpm,cage_adore,cage_harris,heat_W,torque_adore_Nm,torque_palmgren_Nm,torque_skf_Nm")
    for idx in 1:n_speeds
        @printf(f, "%d,%.2f,%.2f,%.4f,%.6f,%.6f,%.6f\n",
            speeds_rpm[idx], cage_speeds_rpm[idx], harris_cage_rpm[idx],
            total_heat_W[idx], friction_torque_Nm[idx], palmgren_torque[idx], skf_torque[idx])
    end
end
println("\n[OK] sweep_data.csv saved")
println("All parametric sweep plots saved to paper/figures/")
