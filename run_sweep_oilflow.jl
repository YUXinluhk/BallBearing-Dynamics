#!/usr/bin/env julia
"""
run_sweep_oilflow.jl — Oil Flow Rate Sensitivity Analysis
Target: Applied Sciences parametric study (Section 4.1)
Bearing: NASA 35-mm bore @ 50,000 RPM (DN = 1.75×10⁶)
Sweep:  V̇ = [50, 100, 200, 400, 760, 1500] cm³/min
Output: Steady-state temperatures, total heat generation
"""

println("="^60)
println("  ADORE — Oil Flow Rate Sensitivity Sweep")
println("="^60)
flush(stdout)

using ADORE
using Printf
using Statistics: mean
using Plots;
ENV["GKSwstype"] = "100"
gr()

default(
    fontfamily="Computer Modern",
    linewidth=1.5,
    framestyle=:box,
    grid=true, gridalpha=0.3,
    titlefontsize=14, guidefontsize=12, tickfontsize=10, legendfontsize=10,
    dpi=300
)

# ─── Base Configuration: NASA 35-mm ───
base_config_file = joinpath(@__DIR__, "inputs", "NASA_35mm_72k_accel.toml")
geom, mat, lub, trac, cage, base_config = load_simulation_config(base_config_file)

# ─── Sweep Parameters ───
oil_flow_rates = [50.0, 100.0, 200.0, 400.0, 760.0, 1500.0]  # cm³/min
target_RPM = 50000.0  # DN = 50000 × 35 = 1.75×10⁶
target_omega = target_RPM * π / 30.0

out_dir = joinpath(@__DIR__, "results", "sweeps", "oilflow")
mkpath(out_dir)

# ─── Storage ───
n_sweep = length(oil_flow_rates)
T_ir_ss  = zeros(n_sweep)
T_or_ss  = zeros(n_sweep)
T_ball_ss = zeros(n_sweep)
T_oil_ss = zeros(n_sweep)
H_total_ss = zeros(n_sweep)

for (idx, V_dot) in enumerate(oil_flow_rates)
    @printf("\n─── Case %d/%d: V̇ = %.0f cm³/min ───\n", idx, n_sweep, V_dot)
    flush(stdout)

    # Rebuild thermal params with new oil flow rate
    th = base_config.thermal
    new_thermal = ThermalParams(
        T_ambient=th.T_ambient,
        T_init=th.T_init,
        T_ref=th.T_ref,
        C_ir=th.C_ir,
        C_or=th.C_or,
        C_ball=th.C_ball,
        C_oil=th.C_oil,
        G_ir_ball=th.G_ir_ball,
        G_or_ball=th.G_or_ball,
        G_ball_oil=th.G_ball_oil,
        G_or_amb=th.G_or_amb,
        G_oil_amb=th.G_oil_amb,
        oil_flow_rate=V_dot,
        T_oil_inlet=th.T_oil_inlet,
    )

    config = SimulationConfig(
        t_end=base_config.t_end,
        dt_output=base_config.dt_output,
        inner_race_speed=target_omega,
        F_axial=base_config.F_axial,
        F_radial=base_config.F_radial,
        t_ramp_end=base_config.t_ramp_end,
        mu_spin=base_config.mu_spin,
        c_structural=base_config.c_structural,
        zeta=base_config.zeta,
        delta_r_thermal=base_config.delta_r_thermal,
        integrator=base_config.integrator,
        churning=base_config.churning,
        thermal=new_thermal,
    )

    result = run_simulation(geom, mat, lub, trac, cage, config; verbose=false)
    @printf("  retcode: %s, steps: %d\n", result.retcode, length(result.t))

    # Extract steady-state temperatures (last 20%)
    Z = geom.n_balls
    u_mat = result.u
    n_t = length(result.t)
    th_off = ADORE.thermal_offset(Z)

    n_ss = max(1, n_t ÷ 5)
    ss_range = max(1, n_t - n_ss + 1):n_t

    T_ir_ss[idx]   = mean(u_mat[ss_range, th_off])
    T_or_ss[idx]   = mean(u_mat[ss_range, th_off+1])
    T_ball_ss[idx] = mean(u_mat[ss_range, th_off+2])
    T_oil_ss[idx]  = mean(u_mat[ss_range, th_off+3])

    # Heat generation from field outputs
    fo_matrix = compute_field_outputs(result)
    fo_ball = zeros(n_t, Z, ADORE.N_FIELD_PER_BALL)
    for j in 1:Z
        base = (j - 1) * ADORE.N_FIELD_PER_BALL
        for f in 1:ADORE.N_FIELD_PER_BALL
            fo_ball[:, j, f] .= fo_matrix[:, base+f]
        end
    end
    H_slide_i = fo_ball[:, :, ADORE.FO_H_SLIDE_I]
    H_slide_o = fo_ball[:, :, ADORE.FO_H_SLIDE_O]
    H_spin_i  = fo_ball[:, :, ADORE.FO_H_SPIN_I]
    H_spin_o  = fo_ball[:, :, ADORE.FO_H_SPIN_O]
    H_drag    = fo_ball[:, :, ADORE.FO_H_DRAG]
    H_churn   = fo_ball[:, :, ADORE.FO_H_CHURN]
    H_total   = H_slide_i .+ H_slide_o .+ H_spin_i .+ H_spin_o .+ H_drag .+ H_churn
    H_bearing = vec(sum(H_total, dims=2))
    H_total_ss[idx] = mean(H_bearing[ss_range])

    @printf("  T_IR=%.1f K, T_OR=%.1f K, T_ball=%.1f K, T_oil=%.1f K\n",
        T_ir_ss[idx], T_or_ss[idx], T_ball_ss[idx], T_oil_ss[idx])
    @printf("  H_total = %.1f W\n", H_total_ss[idx])
    flush(stdout)
end

# ─── Save CSV ───
csv_path = joinpath(out_dir, "oilflow_sweep_results.csv")
open(csv_path, "w") do io
    println(io, "V_dot_cm3min,T_IR_K,T_OR_K,T_ball_K,T_oil_K,H_total_W")
    for i in 1:n_sweep
        @printf(io, "%.1f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            oil_flow_rates[i], T_ir_ss[i], T_or_ss[i], T_ball_ss[i], T_oil_ss[i], H_total_ss[i])
    end
end
println("\n  [OK] CSV saved: $csv_path")

# ─── Publication Figure ───
# Convert to °C
T_ir_C = T_ir_ss .- 273.15
T_or_C = T_or_ss .- 273.15
T_ball_C = T_ball_ss .- 273.15
T_oil_C = T_oil_ss .- 273.15

fig = plot(size=(900, 600), legend=:topright,
    xlabel="Oil Flow Rate (cm³/min)", ylabel="Steady-State Temperature (°C)",
    title="Effect of Oil Flow Rate on Bearing Temperatures\n(NASA 35-mm, 50,000 RPM, Fₐ = 667 N)")

plot!(fig, oil_flow_rates, T_ball_C, marker=:circle, ms=6, color=:darkred, lw=2, label="Ball")
plot!(fig, oil_flow_rates, T_ir_C, marker=:square, ms=5, color=:steelblue, lw=2, label="Inner Race")
plot!(fig, oil_flow_rates, T_or_C, marker=:diamond, ms=5, color=:orange, lw=2, label="Outer Race")
plot!(fig, oil_flow_rates, T_oil_C, marker=:utriangle, ms=5, color=:darkgreen, lw=2, linestyle=:dash, label="Oil Sump")

# Mark NASA baseline
vline!(fig, [760.0], color=:gray, linestyle=:dot, lw=1.5, label="NASA Baseline (760)")

savefig(fig, joinpath(out_dir, "oilflow_temperature_sensitivity.png"))
println("  [OK] Figure saved: oilflow_temperature_sensitivity.png")

# ─── Heat vs Oil Flow ───
fig2 = plot(size=(900, 500), legend=:topright,
    xlabel="Oil Flow Rate (cm³/min)", ylabel="Total Heat Generation (W)",
    title="Heat Generation vs Oil Flow Rate\n(NASA 35-mm, 50,000 RPM)")
plot!(fig2, oil_flow_rates, H_total_ss, marker=:circle, ms=6, color=:black, lw=2, label="H_total")
vline!(fig2, [760.0], color=:gray, linestyle=:dot, lw=1.5, label="NASA Baseline (760)")

savefig(fig2, joinpath(out_dir, "oilflow_heat_generation.png"))
println("  [OK] Figure saved: oilflow_heat_generation.png")

println("\n" * "="^60)
println("  OIL FLOW SWEEP COMPLETE → $out_dir")
println("="^60)
