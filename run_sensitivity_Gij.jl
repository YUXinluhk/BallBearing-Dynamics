#!/usr/bin/env julia
"""
run_sensitivity_Gij.jl — Thermal Conductance Sensitivity Analysis
Perturbs all G_ij values by ±50% to quantify temperature uncertainty
from thermal network parameters (R2 Major Comment M2).
"""

println("="^60)
println("  ADORE — G_ij Sensitivity Analysis (±50%)")
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
    titlefontsize=12, guidefontsize=11, tickfontsize=10, legendfontsize=9,
    dpi=300
)

# ─── Base Configuration ───
base_config_file = joinpath(@__DIR__, "inputs", "NASA_35mm_72k_accel.toml")
geom, mat, lub, trac, cage, base_config = load_simulation_config(base_config_file)

target_RPM = 50000.0
target_omega = target_RPM * π / 30.0

out_dir = joinpath(@__DIR__, "results", "sweeps", "Gij_sensitivity")
mkpath(out_dir)

# ─── Sensitivity Cases ───
# Scale factors for all G values simultaneously
G_scales = [0.5, 0.75, 1.0, 1.25, 1.5]
labels = ["-50%", "-25%", "Baseline", "+25%", "+50%"]

n_cases = length(G_scales)
T_ir_ss  = zeros(n_cases)
T_or_ss  = zeros(n_cases)
T_ball_ss = zeros(n_cases)
T_oil_ss = zeros(n_cases)
H_total_ss = zeros(n_cases)
dT_ir_or = zeros(n_cases)

for (idx, scale) in enumerate(G_scales)
    @printf("\n─── Case %d/%d: G × %.2f (%s) ───\n", idx, n_cases, scale, labels[idx])
    flush(stdout)

    th = base_config.thermal
    new_thermal = ThermalParams(
        T_ambient=th.T_ambient,
        T_init=th.T_init,
        T_ref=th.T_ref,
        C_ir=th.C_ir,
        C_or=th.C_or,
        C_ball=th.C_ball,
        C_oil=th.C_oil,
        G_ir_ball=th.G_ir_ball * scale,
        G_or_ball=th.G_or_ball * scale,
        G_ball_oil=th.G_ball_oil * scale,
        G_or_amb=th.G_or_amb * scale,
        G_oil_amb=th.G_oil_amb * scale,
        oil_flow_rate=th.oil_flow_rate,
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

    Z = geom.n_balls
    u_mat = result.u
    n_t = length(result.t)
    th_off = ADORE.thermal_offset(Z)
    n_ss = max(1, n_t ÷ 5)
    ss_range = max(1, n_t - n_ss + 1):n_t

    T_ir_ss[idx]  = mean(u_mat[ss_range, th_off])
    T_or_ss[idx]  = mean(u_mat[ss_range, th_off+1])
    T_ball_ss[idx] = mean(u_mat[ss_range, th_off+2])
    T_oil_ss[idx] = mean(u_mat[ss_range, th_off+3])
    dT_ir_or[idx] = T_ir_ss[idx] - T_or_ss[idx]

    fo_matrix = compute_field_outputs(result)
    fo_ball = zeros(n_t, Z, ADORE.N_FIELD_PER_BALL)
    for j in 1:Z
        base = (j - 1) * ADORE.N_FIELD_PER_BALL
        for f in 1:ADORE.N_FIELD_PER_BALL
            fo_ball[:, j, f] .= fo_matrix[:, base+f]
        end
    end
    H_total = vec(sum(
        fo_ball[:, :, ADORE.FO_H_SLIDE_I] .+
        fo_ball[:, :, ADORE.FO_H_SLIDE_O] .+
        fo_ball[:, :, ADORE.FO_H_SPIN_I] .+
        fo_ball[:, :, ADORE.FO_H_SPIN_O] .+
        fo_ball[:, :, ADORE.FO_H_DRAG] .+
        fo_ball[:, :, ADORE.FO_H_CHURN], dims=2))
    H_total_ss[idx] = mean(H_total[ss_range])

    @printf("  T_IR=%.2f K, T_OR=%.2f K, ΔT=%.2f K, H=%.1f W\n",
        T_ir_ss[idx], T_or_ss[idx], dT_ir_or[idx], H_total_ss[idx])
    flush(stdout)
end

# ─── Save CSV ───
csv_path = joinpath(out_dir, "Gij_sensitivity_results.csv")
open(csv_path, "w") do io
    println(io, "G_scale,label,T_IR_K,T_OR_K,T_ball_K,T_oil_K,dT_IR_OR_K,H_total_W")
    for i in 1:n_cases
        @printf(io, "%.2f,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            G_scales[i], labels[i], T_ir_ss[i], T_or_ss[i], T_ball_ss[i],
            T_oil_ss[i], dT_ir_or[i], H_total_ss[i])
    end
end
println("\n  [OK] CSV saved: $csv_path")

# ─── Bar chart figure ───
fig = plot(layout=(1,2), size=(1200, 500), margin=5Plots.mm)

# Panel a: Temperatures
x = 1:n_cases
plot!(fig[1], x, T_ir_ss .- 273.15, marker=:circle, ms=6, color=:steelblue, lw=2, label="T_IR",
    xticks=(x, labels), xlabel="G_ij perturbation", ylabel="Temperature (°C)",
    title="(a) Steady-state temperatures")
plot!(fig[1], x, T_or_ss .- 273.15, marker=:square, ms=5, color=:orange, lw=2, label="T_OR")
plot!(fig[1], x, T_ball_ss .- 273.15, marker=:diamond, ms=5, color=:darkred, lw=2, label="T_ball")

# Panel b: ΔT_IR-OR
plot!(fig[2], x, dT_ir_or, marker=:circle, ms=6, color=:black, lw=2, label="ΔT_IR-OR",
    xticks=(x, labels), xlabel="G_ij perturbation", ylabel="ΔT (K)",
    title="(b) Inner-outer race differential")

savefig(fig, joinpath(out_dir, "Gij_sensitivity.png"))
println("  [OK] Figure saved: Gij_sensitivity.png")

# ─── Summary ───
base_idx = 3  # baseline
println("\n" * "="^60)
println("  SENSITIVITY SUMMARY (relative to baseline)")
println("="^60)
for i in 1:n_cases
    i == base_idx && continue
    @printf("  %s: ΔT_IR = %+.2f K, ΔH = %+.1f W\n",
        labels[i], T_ir_ss[i] - T_ir_ss[base_idx], H_total_ss[i] - H_total_ss[base_idx])
end
println("="^60)
