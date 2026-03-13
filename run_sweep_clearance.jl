#!/usr/bin/env julia
"""
run_sweep_clearance.jl — Diametral Clearance vs Thermal Expansion Analysis
Target: Applied Sciences parametric study (Section 4.2)
Bearing: NASA 35-mm bore @ 50,000 RPM
Sweep:  P_d = [0, 10, 20, 40, 60, 100] μm
Loading: F_axial = 100 N (light preload), F_radial = 500 N
Focus:  Show how clearance affects contact angle and temperature
        under combined axial/radial loading
"""

println("="^60)
println("  ADORE — Clearance Sensitivity Sweep")
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
geom_base, mat, lub, trac, cage_base, base_config = load_simulation_config(base_config_file)

# ─── Sweep Parameters ───
clearances_um = [0.0, 10.0, 20.0, 40.0, 60.0, 100.0]  # μm
clearances_m = clearances_um .* 1e-6
target_RPM = 50000.0
target_omega = target_RPM * π / 30.0

out_dir = joinpath(@__DIR__, "results", "sweeps", "clearance")
mkpath(out_dir)

n_sweep = length(clearances_um)

# Store full transient data for overlay plot
all_t_ms = Vector{Vector{Float64}}(undef, n_sweep)
all_T_ir = Vector{Vector{Float64}}(undef, n_sweep)
all_T_or = Vector{Vector{Float64}}(undef, n_sweep)
all_T_ball = Vector{Vector{Float64}}(undef, n_sweep)
all_alpha_i_mean = Vector{Vector{Float64}}(undef, n_sweep)  # mean inner contact angle
T_ir_ss = zeros(n_sweep)
T_or_ss = zeros(n_sweep)
T_ball_ss = zeros(n_sweep)
H_total_ss = zeros(n_sweep)

for (idx, Pd) in enumerate(clearances_m)
    @printf("\n─── Case %d/%d: P_d = %.1f μm ───\n", idx, n_sweep, clearances_um[idx])
    flush(stdout)

    # Rebuild geometry with new clearance
    geom = BearingGeometry(
        d=geom_base.d,
        n_balls=geom_base.n_balls,
        f_i=geom_base.f_i,
        f_o=geom_base.f_o,
        d_m=geom_base.d_m,
        alpha_0=geom_base.alpha_0,
        P_d=Pd,
        rho_ball=geom_base.rho_ball,
    )
    cage = cage_from_bearing(geom;
        pocket_clearance_ratio=0.01,
        pilot_clearance_ratio=0.002)

    config = SimulationConfig(
        t_end=base_config.t_end,
        dt_output=base_config.dt_output,
        inner_race_speed=target_omega,
        F_axial=100.0,       # Light axial preload (down from 667 N)
        F_radial=500.0,      # Radial load → clearance-sensitive load zone
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

    # Full transient
    all_t_ms[idx] = result.t .* 1e3
    all_T_ir[idx] = u_mat[:, th_off] .- 273.15
    all_T_or[idx] = u_mat[:, th_off+1] .- 273.15
    all_T_ball[idx] = u_mat[:, th_off+2] .- 273.15

    # Compute field outputs (contact angles, loads, heat gen)
    fo_matrix = compute_field_outputs(result)
    fo_ball = zeros(n_t, Z, ADORE.N_FIELD_PER_BALL)
    for j in 1:Z
        base_fo = (j - 1) * ADORE.N_FIELD_PER_BALL
        for f in 1:ADORE.N_FIELD_PER_BALL
            fo_ball[:, j, f] .= fo_matrix[:, base_fo+f]
        end
    end

    # Mean inner contact angle (averaged over all balls), in degrees
    alpha_i_all = fo_ball[:, :, ADORE.FO_ALPHA_I]  # n_t × Z matrix
    all_alpha_i_mean[idx] = vec(mean(alpha_i_all, dims=2)) .* (180.0 / π)

    # Steady state
    n_ss = max(1, n_t ÷ 5)
    ss_range = max(1, n_t - n_ss + 1):n_t
    T_ir_ss[idx]   = mean(u_mat[ss_range, th_off])
    T_or_ss[idx]   = mean(u_mat[ss_range, th_off+1])
    T_ball_ss[idx] = mean(u_mat[ss_range, th_off+2])

    # Heat generation
    H_total = fo_ball[:, :, ADORE.FO_H_SLIDE_I] .+ fo_ball[:, :, ADORE.FO_H_SLIDE_O] .+
              fo_ball[:, :, ADORE.FO_H_SPIN_I] .+ fo_ball[:, :, ADORE.FO_H_SPIN_O] .+
              fo_ball[:, :, ADORE.FO_H_DRAG] .+ fo_ball[:, :, ADORE.FO_H_CHURN]
    H_bearing = vec(sum(H_total, dims=2))
    H_total_ss[idx] = mean(H_bearing[ss_range])

    @printf("  T_IR=%.1f K, T_ball=%.1f K, α_i=%.2f°, H=%.1f W\n",
        T_ir_ss[idx], T_ball_ss[idx], all_alpha_i_mean[idx][end], H_total_ss[idx])
    flush(stdout)
end

# ─── Save CSV ───
csv_path = joinpath(out_dir, "clearance_sweep_results.csv")
open(csv_path, "w") do io
    println(io, "P_d_um,T_IR_K,T_OR_K,T_ball_K,H_total_W")
    for i in 1:n_sweep
        @printf(io, "%.1f,%.3f,%.3f,%.3f,%.3f\n",
            clearances_um[i], T_ir_ss[i], T_or_ss[i], T_ball_ss[i], H_total_ss[i])
    end
end
println("\n  [OK] CSV: $csv_path")

# ─── Save Transient CSVs (interpolated to common time grid) ───
# Find the shortest time span across all cases, then create a uniform grid
t_max = minimum(last.(all_t_ms))
t_min = maximum(first.(all_t_ms))
n_points = 500  # enough for a smooth plot
t_common = range(t_min, t_max, length=n_points) |> collect

# Helper: linear interpolation
function interp1(t_data, y_data, t_val)
    j = searchsortedlast(t_data, t_val)
    j = clamp(j, 1, length(t_data) - 1)
    frac = (t_val - t_data[j]) / (t_data[j+1] - t_data[j] + 1e-30)
    frac = clamp(frac, 0.0, 1.0)
    return y_data[j] + frac * (y_data[j+1] - y_data[j])
end

# --- Temperature transient CSV ---
transient_csv = joinpath(out_dir, "clearance_transient_data.csv")
open(transient_csv, "w") do io
    header = "time_ms"
    for i in 1:n_sweep
        header *= @sprintf(",T_IR_Pd%.0fum", clearances_um[i])
    end
    println(io, header)
    for k in 1:n_points
        t_val = t_common[k]
        @printf(io, "%.4f", t_val)
        for i in 1:n_sweep
            @printf(io, ",%.4f", interp1(all_t_ms[i], all_T_ir[i], t_val))
        end
        println(io)
    end
end
println("  [OK] Transient CSV (T_IR): $transient_csv")

# --- Contact angle transient CSV ---
alpha_csv = joinpath(out_dir, "clearance_alpha_transient_data.csv")
open(alpha_csv, "w") do io
    header = "time_ms"
    for i in 1:n_sweep
        header *= @sprintf(",alpha_i_Pd%.0fum", clearances_um[i])
    end
    println(io, header)
    for k in 1:n_points
        t_val = t_common[k]
        @printf(io, "%.4f", t_val)
        for i in 1:n_sweep
            @printf(io, ",%.4f", interp1(all_t_ms[i], all_alpha_i_mean[i], t_val))
        end
        println(io)
    end
end
println("  [OK] Transient CSV (α_i): $alpha_csv")

# ─── Figure 1: Contact angle transient overlay ───
colors_sweep = [:darkred, :red, :orange, :steelblue, :blue, :purple]

fig1 = plot(size=(1000, 600), legend=:bottomright,
    xlabel="Time (ms)", ylabel="Mean Inner Contact Angle (°)",
    title="Effect of Initial Clearance on Contact Angle Transient\n(NASA 35-mm, 50,000 RPM)")

for i in 1:n_sweep
    plot!(fig1, all_t_ms[i], all_alpha_i_mean[i],
        color=colors_sweep[i], lw=1.8,
        label=@sprintf("Pd = %.0f μm", clearances_um[i]))
end

savefig(fig1, joinpath(out_dir, "clearance_transient_overlay.png"))
println("  [OK] Figure: clearance_transient_overlay.png")

# ─── Figure 2: Steady-state T_ball and H vs clearance ───
fig2 = plot(layout=(1, 2), size=(1200, 500))

plot!(fig2[1], clearances_um, T_ball_ss .- 273.15,
    marker=:circle, ms=6, color=:darkred, lw=2, label="Ball",
    xlabel="Diametral Clearance (μm)", ylabel="Steady-State Temperature (°C)",
    title="Temperature vs Clearance")
plot!(fig2[1], clearances_um, T_ir_ss .- 273.15,
    marker=:square, ms=5, color=:steelblue, lw=2, label="Inner Race")
plot!(fig2[1], clearances_um, T_or_ss .- 273.15,
    marker=:diamond, ms=5, color=:orange, lw=2, label="Outer Race")

plot!(fig2[2], clearances_um, H_total_ss,
    marker=:circle, ms=6, color=:black, lw=2, label="",
    xlabel="Diametral Clearance (μm)", ylabel="Total Heat Generation (W)",
    title="Heat Generation vs Clearance")

savefig(fig2, joinpath(out_dir, "clearance_steadystate.png"))
println("  [OK] Figure: clearance_steadystate.png")

println("\n" * "="^60)
println("  CLEARANCE SWEEP COMPLETE → $out_dir")
println("="^60)
