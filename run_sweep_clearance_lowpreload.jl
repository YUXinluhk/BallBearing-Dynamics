#!/usr/bin/env julia
"""
run_sweep_clearance_lowpreload.jl — Clearance Sweep Under Low Axial Preload
Bearing: NASA 35-mm bore @ 50,000 RPM, F_a = 100 N (light preload)
Sweep:  P_d = [0, 10, 20, 30, 50, 80, 100] μm
Addresses R2 Major Comment M4: show clearance sensitivity under different loading.
"""

println("="^60)
println("  ADORE — Clearance Sweep (Low Preload, F_a = 100 N)")
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
F_axial_low = 100.0  # Low preload

# ─── Sweep Parameters ───
clearance_um = [0.0, 10.0, 20.0, 30.0, 50.0, 80.0, 100.0]  # μm
clearance_m = clearance_um .* 1e-6  # Convert to meters

out_dir = joinpath(@__DIR__, "results", "sweeps", "clearance_lowpreload")
mkpath(out_dir)

n_sweep = length(clearance_um)
T_ir_ss  = zeros(n_sweep)
T_or_ss  = zeros(n_sweep)
T_ball_ss = zeros(n_sweep)
T_oil_ss = zeros(n_sweep)
H_total_ss = zeros(n_sweep)

for (idx, Pd) in enumerate(clearance_m)
    @printf("\n─── Case %d/%d: P_d = %.0f μm, F_a = %.0f N ───\n",
        idx, n_sweep, clearance_um[idx], F_axial_low)
    flush(stdout)

    # Rebuild geometry with new clearance
    new_geom = BearingGeometry(
        d=geom.d,
        n_balls=geom.n_balls,
        f_i=geom.f_i,
        f_o=geom.f_o,
        d_m=geom.d_m,
        alpha_0=geom.alpha_0,
        P_d=Pd,
        rho_ball=geom.rho_ball,
    )

    config = SimulationConfig(
        t_end=base_config.t_end,
        dt_output=base_config.dt_output,
        inner_race_speed=target_omega,
        F_axial=F_axial_low,
        F_radial=base_config.F_radial,
        t_ramp_end=base_config.t_ramp_end,
        mu_spin=base_config.mu_spin,
        c_structural=base_config.c_structural,
        zeta=base_config.zeta,
        delta_r_thermal=base_config.delta_r_thermal,
        integrator=base_config.integrator,
        churning=base_config.churning,
        thermal=base_config.thermal,
    )

    result = run_simulation(new_geom, mat, lub, trac, cage, config; verbose=false)
    @printf("  retcode: %s, steps: %d\n", result.retcode, length(result.t))

    Z = new_geom.n_balls
    u_mat = result.u
    n_t = length(result.t)
    th_off = ADORE.thermal_offset(Z)
    n_ss = max(1, n_t ÷ 5)
    ss_range = max(1, n_t - n_ss + 1):n_t

    T_ir_ss[idx]  = mean(u_mat[ss_range, th_off])
    T_or_ss[idx]  = mean(u_mat[ss_range, th_off+1])
    T_ball_ss[idx] = mean(u_mat[ss_range, th_off+2])
    T_oil_ss[idx] = mean(u_mat[ss_range, th_off+3])

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

    @printf("  T_IR=%.2f K, T_OR=%.2f K, H=%.1f W\n",
        T_ir_ss[idx], T_or_ss[idx], H_total_ss[idx])
    flush(stdout)
end

# ─── Save CSV ───
csv_path = joinpath(out_dir, "clearance_lowpreload_results.csv")
open(csv_path, "w") do io
    println(io, "Pd_um,T_IR_K,T_OR_K,T_ball_K,T_oil_K,H_total_W")
    for i in 1:n_sweep
        @printf(io, "%.0f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            clearance_um[i], T_ir_ss[i], T_or_ss[i], T_ball_ss[i], T_oil_ss[i], H_total_ss[i])
    end
end
println("\n  [OK] CSV saved: $csv_path")

# ─── Publication Figure ───
fig = plot(layout=(1,2), size=(1200, 500), margin=5Plots.mm)

plot!(fig[1], clearance_um, T_ir_ss .- 273.15,
    marker=:circle, ms=6, color=:steelblue, lw=2, label="T_IR",
    xlabel="Initial Clearance P_d (μm)", ylabel="Temperature (°C)",
    title="(a) Low preload (F_a = 100 N)")
plot!(fig[1], clearance_um, T_or_ss .- 273.15,
    marker=:square, ms=5, color=:orange, lw=2, label="T_OR")
plot!(fig[1], clearance_um, T_ball_ss .- 273.15,
    marker=:diamond, ms=5, color=:darkred, lw=2, label="T_ball")

plot!(fig[2], clearance_um, H_total_ss,
    marker=:circle, ms=6, color=:black, lw=2, label="H_total",
    xlabel="Initial Clearance P_d (μm)", ylabel="Heat Generation (W)",
    title="(b) Total heat generation")

savefig(fig, joinpath(out_dir, "clearance_lowpreload.png"))
println("  [OK] Figure saved: clearance_lowpreload.png")

println("\n" * "="^60)
println("  CLEARANCE SWEEP (LOW PRELOAD) COMPLETE → $out_dir")
println("="^60)
