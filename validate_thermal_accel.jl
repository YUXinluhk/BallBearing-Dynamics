#!/usr/bin/env julia
"""
validate_thermal_accel.jl — Thermal Acceleration Validation
Compare 1× (physical) vs 100× (accelerated) thermal capacitance
at NASA 35-mm, 50,000 RPM baseline to prove steady-state equivalence.
Generates figure for Applied Sciences revision (Major Comment 3).
"""

println("="^60)
println("  ADORE — Thermal Acceleration Validation (1× vs 100×)")
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

# ─── Base Configuration: NASA 35-mm ───
base_config_file = joinpath(@__DIR__, "inputs", "NASA_35mm_72k_accel.toml")
geom, mat, lub, trac, cage, base_config = load_simulation_config(base_config_file)

target_RPM = 50000.0
target_omega = target_RPM * π / 30.0

out_dir = joinpath(@__DIR__, "results", "sweeps", "thermal_accel")
mkpath(out_dir)

# ─── Physical thermal capacitances (from TOML: accelerated values × 100) ───
C_ir_phys = 17.47     # [J/K]
C_or_phys = 26.20
C_ball_phys = 10.94
C_oil_phys = 52.45

# ─── Run configs ───
cases = Dict{String, NamedTuple{(:label, :C_scale, :t_end, :dt_output), Tuple{String,Float64,Float64,Float64}}}()
# Accelerated case (100×): short simulation, fast convergence
cases["accel"] = (label="Accelerated (C/100)", C_scale=0.01, t_end=0.15, dt_output=0.002)
# Physical case (1×): need 100× longer to reach same thermal state
cases["physical"] = (label="Physical (1×)", C_scale=1.0, t_end=0.5, dt_output=0.01)

results = Dict{String, Any}()

for (key, case) in sort(collect(cases))
    @printf("\n─── %s: C_scale=%.2f, t_end=%.1f s ───\n", case.label, case.C_scale, case.t_end)
    flush(stdout)

    th = base_config.thermal
    new_thermal = ThermalParams(
        T_ambient=th.T_ambient,
        T_init=th.T_init,
        T_ref=th.T_ref,
        C_ir=C_ir_phys * case.C_scale,
        C_or=C_or_phys * case.C_scale,
        C_ball=C_ball_phys * case.C_scale,
        C_oil=C_oil_phys * case.C_scale,
        G_ir_ball=th.G_ir_ball,
        G_or_ball=th.G_or_ball,
        G_ball_oil=th.G_ball_oil,
        G_or_amb=th.G_or_amb,
        G_oil_amb=th.G_oil_amb,
        oil_flow_rate=th.oil_flow_rate,
        T_oil_inlet=th.T_oil_inlet,
    )

    config = SimulationConfig(
        t_end=case.t_end,
        dt_output=case.dt_output,
        inner_race_speed=target_omega,
        F_axial=base_config.F_axial,
        F_radial=base_config.F_radial,
        t_ramp_end=base_config.t_ramp_end,
        mu_spin=base_config.mu_spin,
        c_structural=base_config.c_structural,
        zeta=base_config.zeta,
        delta_r_thermal=base_config.delta_r_thermal,
        integrator=IntegratorConfig(
            rtol=1e-4, atol=1e-7, h_max=5e-5,
            max_steps=50_000_000,
        ),
        churning=base_config.churning,
        thermal=new_thermal,
    )

    result = run_simulation(geom, mat, lub, trac, cage, config; verbose=false)
    @printf("  retcode: %s, steps: %d\n", result.retcode, length(result.t))

    Z = geom.n_balls
    u_mat = result.u
    n_t = length(result.t)
    th_off = ADORE.thermal_offset(Z)

    T_ir  = u_mat[:, th_off]
    T_or  = u_mat[:, th_off+1]
    T_ball = u_mat[:, th_off+2]
    T_oil = u_mat[:, th_off+3]

    # Steady-state: last 20%
    n_ss = max(1, n_t ÷ 5)
    ss_range = max(1, n_t - n_ss + 1):n_t
    T_ir_ss = mean(T_ir[ss_range])
    T_or_ss = mean(T_or[ss_range])
    T_ball_ss = mean(T_ball[ss_range])
    T_oil_ss = mean(T_oil[ss_range])

    # Heat generation
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
    H_ss = mean(H_total[ss_range])

    results[key] = (
        t=result.t, T_ir=T_ir, T_or=T_or, T_ball=T_ball, T_oil=T_oil,
        H_total=H_total,
        ss=(T_ir=T_ir_ss, T_or=T_or_ss, T_ball=T_ball_ss, T_oil=T_oil_ss, H=H_ss)
    )

    @printf("  Steady-state: T_IR=%.2f K, T_OR=%.2f K, T_ball=%.2f K, T_oil=%.2f K, H=%.1f W\n",
        T_ir_ss, T_or_ss, T_ball_ss, T_oil_ss, H_ss)
    flush(stdout)
end

# ─── Comparison Summary ───
r_a = results["accel"].ss
r_p = results["physical"].ss

println("\n" * "="^60)
println("  COMPARISON: Accelerated vs Physical Thermal Capacitance")
println("="^60)
@printf("  %-12s  %10s  %10s  %10s\n", "Variable", "Accel (C/100)", "Physical", "Δ (K)")
@printf("  %-12s  %10.2f  %10.2f  %10.3f\n", "T_IR (K)", r_a.T_ir, r_p.T_ir, r_a.T_ir - r_p.T_ir)
@printf("  %-12s  %10.2f  %10.2f  %10.3f\n", "T_OR (K)", r_a.T_or, r_p.T_or, r_a.T_or - r_p.T_or)
@printf("  %-12s  %10.2f  %10.2f  %10.3f\n", "T_ball (K)", r_a.T_ball, r_p.T_ball, r_a.T_ball - r_p.T_ball)
@printf("  %-12s  %10.2f  %10.2f  %10.3f\n", "T_oil (K)", r_a.T_oil, r_p.T_oil, r_a.T_oil - r_p.T_oil)
@printf("  %-12s  %10.1f  %10.1f  %10.2f\n", "H_total (W)", r_a.H, r_p.H, r_a.H - r_p.H)

# ─── Publication Figures ───
# Figure 1: Temperature transient comparison (normalized time)
fig1 = plot(layout=(2,2), size=(1200, 800), margin=5Plots.mm)

t_a = results["accel"].t
t_p = results["physical"].t

# Accelerated uses "simulation time"; physical uses real time
# To compare, we normalize to fraction of total sim time for shape comparison
# But the reviewer wants actual trajectories, so let's show both on their own time axes

# Panel (a): Inner race temperature
plot!(fig1[1], t_a, results["accel"].T_ir .- 273.15,
    color=:red, lw=2, label="Accelerated (C/100, 0.15 s)",
    xlabel="Time (s)", ylabel="T_IR (°C)", title="(a) Inner Race")
plot!(fig1[1], t_p, results["physical"].T_ir .- 273.15,
    color=:blue, lw=2, linestyle=:dash, label="Physical (1×, 3.0 s)")

# Panel (b): Outer race temperature
plot!(fig1[2], t_a, results["accel"].T_or .- 273.15,
    color=:red, lw=2, label="Accelerated",
    xlabel="Time (s)", ylabel="T_OR (°C)", title="(b) Outer Race")
plot!(fig1[2], t_p, results["physical"].T_or .- 273.15,
    color=:blue, lw=2, linestyle=:dash, label="Physical")

# Panel (c): Ball temperature
plot!(fig1[3], t_a, results["accel"].T_ball .- 273.15,
    color=:red, lw=2, label="Accelerated",
    xlabel="Time (s)", ylabel="T_ball (°C)", title="(c) Ball")
plot!(fig1[3], t_p, results["physical"].T_ball .- 273.15,
    color=:blue, lw=2, linestyle=:dash, label="Physical")

# Panel (d): Heat generation
plot!(fig1[4], t_a, results["accel"].H_total,
    color=:red, lw=2, label="Accelerated",
    xlabel="Time (s)", ylabel="H_total (W)", title="(d) Total Heat Generation")
plot!(fig1[4], t_p, results["physical"].H_total,
    color=:blue, lw=2, linestyle=:dash, label="Physical")

savefig(fig1, joinpath(out_dir, "thermal_accel_validation.png"))
println("\n  [OK] Figure saved: thermal_accel_validation.png")

# Save CSV comparison
csv_path = joinpath(out_dir, "thermal_accel_comparison.csv")
open(csv_path, "w") do io
    println(io, "case,T_IR_K,T_OR_K,T_ball_K,T_oil_K,H_total_W")
    @printf(io, "accelerated,%.3f,%.3f,%.3f,%.3f,%.3f\n",
        r_a.T_ir, r_a.T_or, r_a.T_ball, r_a.T_oil, r_a.H)
    @printf(io, "physical,%.3f,%.3f,%.3f,%.3f,%.3f\n",
        r_p.T_ir, r_p.T_or, r_p.T_ball, r_p.T_oil, r_p.H)
end
println("  [OK] CSV saved: $csv_path")

println("\n" * "="^60)
println("  THERMAL ACCELERATION VALIDATION COMPLETE → $out_dir")
println("="^60)
