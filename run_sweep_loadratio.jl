#!/usr/bin/env julia
"""
run_sweep_loadratio.jl — Axial/Radial Load Ratio vs Temperature Split
Target: Applied Sciences parametric study (Section 4.3)
Bearing: NASA 35-mm bore @ 50,000 RPM
Sweep:  Fa/Fr ratio from pure-radial to pure-axial (constant |P| = 1000 N)
Output: IR-OR temperature difference, heat source redistribution
"""

println("="^60)
println("  ADORE — Load Ratio Sensitivity Sweep")
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

# ─── Base Configuration ───
base_config_file = joinpath(@__DIR__, "inputs", "NASA_35mm_72k_accel.toml")
geom, mat, lub, trac, cage, base_config = load_simulation_config(base_config_file)

# ─── Sweep Parameters ───
P_total = 1000.0  # N — constant resultant load magnitude
# Load angles: 0° = pure radial, 90° = pure axial
load_angles_deg = [10.0, 25.0, 45.0, 65.0, 80.0, 90.0]
target_RPM = 50000.0
target_omega = target_RPM * π / 30.0

out_dir = joinpath(@__DIR__, "results", "sweeps", "loadratio")
mkpath(out_dir)

n_sweep = length(load_angles_deg)
Fa_vals = P_total .* sind.(load_angles_deg)
Fr_vals = P_total .* cosd.(load_angles_deg)

# Storage
T_ir_ss = zeros(n_sweep)
T_or_ss = zeros(n_sweep)
T_ball_ss = zeros(n_sweep)
dT_ir_or = zeros(n_sweep)
H_total_ss = zeros(n_sweep)
H_slide_i_ss = zeros(n_sweep)
H_slide_o_ss = zeros(n_sweep)
H_spin_i_ss = zeros(n_sweep)
H_spin_o_ss = zeros(n_sweep)
H_drag_ss = zeros(n_sweep)
H_churn_ss = zeros(n_sweep)

for (idx, angle) in enumerate(load_angles_deg)
    @printf("\n─── Case %d/%d: θ=%.0f° → Fa=%.0f N, Fr=%.0f N ───\n",
        idx, n_sweep, angle, Fa_vals[idx], Fr_vals[idx])
    flush(stdout)

    config = SimulationConfig(
        t_end=base_config.t_end,
        dt_output=base_config.dt_output,
        inner_race_speed=target_omega,
        F_axial=Fa_vals[idx],
        F_radial=Fr_vals[idx],
        t_ramp_end=base_config.t_ramp_end,
        mu_spin=base_config.mu_spin,
        c_structural=base_config.c_structural,
        zeta=base_config.zeta,
        delta_r_thermal=base_config.delta_r_thermal,
        integrator=base_config.integrator,
        churning=base_config.churning,
        thermal=base_config.thermal,
    )

    result = run_simulation(geom, mat, lub, trac, cage, config; verbose=false)
    @printf("  retcode: %s, steps: %d\n", result.retcode, length(result.t))

    Z = geom.n_balls
    u_mat = result.u
    n_t = length(result.t)
    th_off = ADORE.thermal_offset(Z)

    n_ss = max(1, n_t ÷ 5)
    ss_range = max(1, n_t - n_ss + 1):n_t

    T_ir_ss[idx]   = mean(u_mat[ss_range, th_off])
    T_or_ss[idx]   = mean(u_mat[ss_range, th_off+1])
    T_ball_ss[idx] = mean(u_mat[ss_range, th_off+2])
    dT_ir_or[idx]  = T_ir_ss[idx] - T_or_ss[idx]

    # Heat sources breakdown
    fo_matrix = compute_field_outputs(result)
    fo_ball = zeros(n_t, Z, ADORE.N_FIELD_PER_BALL)
    for j in 1:Z
        base_fo = (j - 1) * ADORE.N_FIELD_PER_BALL
        for f in 1:ADORE.N_FIELD_PER_BALL
            fo_ball[:, j, f] .= fo_matrix[:, base_fo+f]
        end
    end

    H_si = fo_ball[:, :, ADORE.FO_H_SLIDE_I]
    H_so = fo_ball[:, :, ADORE.FO_H_SLIDE_O]
    H_pi = fo_ball[:, :, ADORE.FO_H_SPIN_I]
    H_po = fo_ball[:, :, ADORE.FO_H_SPIN_O]
    H_dr = fo_ball[:, :, ADORE.FO_H_DRAG]
    H_ch = fo_ball[:, :, ADORE.FO_H_CHURN]

    H_slide_i_ss[idx] = mean(sum(H_si[ss_range, :], dims=2))
    H_slide_o_ss[idx] = mean(sum(H_so[ss_range, :], dims=2))
    H_spin_i_ss[idx]  = mean(sum(H_pi[ss_range, :], dims=2))
    H_spin_o_ss[idx]  = mean(sum(H_po[ss_range, :], dims=2))
    H_drag_ss[idx]    = mean(sum(H_dr[ss_range, :], dims=2))
    H_churn_ss[idx]   = mean(sum(H_ch[ss_range, :], dims=2))
    H_total_ss[idx]   = H_slide_i_ss[idx] + H_slide_o_ss[idx] +
                         H_spin_i_ss[idx] + H_spin_o_ss[idx] +
                         H_drag_ss[idx] + H_churn_ss[idx]

    @printf("  T_IR=%.1f K, T_OR=%.1f K, ΔT=%.2f K, H=%.1f W\n",
        T_ir_ss[idx], T_or_ss[idx], dT_ir_or[idx], H_total_ss[idx])
    flush(stdout)
end

# ─── Save CSV ───
csv_path = joinpath(out_dir, "loadratio_sweep_results.csv")
open(csv_path, "w") do io
    println(io, "angle_deg,Fa_N,Fr_N,T_IR_K,T_OR_K,T_ball_K,dT_IR_OR_K,H_total_W,H_slide_i_W,H_slide_o_W,H_spin_i_W,H_spin_o_W,H_drag_W,H_churn_W")
    for i in 1:n_sweep
        @printf(io, "%.1f,%.1f,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            load_angles_deg[i], Fa_vals[i], Fr_vals[i],
            T_ir_ss[i], T_or_ss[i], T_ball_ss[i], dT_ir_or[i], H_total_ss[i],
            H_slide_i_ss[i], H_slide_o_ss[i], H_spin_i_ss[i], H_spin_o_ss[i],
            H_drag_ss[i], H_churn_ss[i])
    end
end
println("\n  [OK] CSV: $csv_path")

# ─── Figure 1: Temperature vs Load Angle ───
Fa_Fr_ratio = Fa_vals ./ max.(Fr_vals, 0.01)

fig1 = plot(layout=(1, 2), size=(1400, 550))

plot!(fig1[1], load_angles_deg, T_ir_ss .- 273.15,
    marker=:square, ms=5, color=:steelblue, lw=2, label="Inner Race",
    xlabel="Load Angle θ (°)\n(0° = radial, 90° = axial)",
    ylabel="Steady-State Temperature (°C)",
    title="Temperature vs Load Direction")
plot!(fig1[1], load_angles_deg, T_or_ss .- 273.15,
    marker=:diamond, ms=5, color=:orange, lw=2, label="Outer Race")
plot!(fig1[1], load_angles_deg, T_ball_ss .- 273.15,
    marker=:circle, ms=6, color=:darkred, lw=2, label="Ball")

plot!(fig1[2], load_angles_deg, dT_ir_or,
    marker=:circle, ms=6, color=:purple, lw=2.5, label="",
    xlabel="Load Angle θ (°)", ylabel="ΔT (IR − OR) [K]",
    title="Inner-Outer Race Temperature Difference")

savefig(fig1, joinpath(out_dir, "loadratio_temperature.png"))
println("  [OK] Figure: loadratio_temperature.png")

# ─── Figure 2: Heat source breakdown stacked bar ───
fig2 = plot(size=(1000, 600), legend=:topright,
    xlabel="Load Angle θ (°)", ylabel="Heat Generation (W)",
    title="Heat Source Breakdown vs Load Direction\n(NASA 35-mm, 50,000 RPM, |P| = 1000 N)")

bar_bottom = zeros(n_sweep)
source_data = [
    (H_slide_i_ss, "Inner Sliding", "#1f77b4"),
    (H_slide_o_ss, "Outer Sliding", "#d62728"),
    (H_spin_i_ss, "Inner Spin", "#2ca02c"),
    (H_spin_o_ss, "Outer Spin", "#ff7f0e"),
    (H_drag_ss, "Churning Drag", "#9467bd"),
    (H_churn_ss, "Churning Rot.", "#8c564b"),
]

for (vals, label, color) in source_data
    bar!(fig2, load_angles_deg, vals, fillrange=bar_bottom, label=label,
        color=color, alpha=0.85, bar_width=8)
    bar_bottom .+= vals
end

savefig(fig2, joinpath(out_dir, "loadratio_heat_breakdown.png"))
println("  [OK] Figure: loadratio_heat_breakdown.png")

println("\n" * "="^60)
println("  LOAD RATIO SWEEP COMPLETE → $out_dir")
println("="^60)
