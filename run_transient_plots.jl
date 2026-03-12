#!/usr/bin/env julia
"""
run_transient_plots.jl — Extended transient simulation + publication-quality plots.

Generates 3 transient dynamics figures for manuscript revision:
  15_Spin_Nutation.png   — Ball spin-axis nutation settling trajectory
  16_Cage_Whirl.png      — Cage center-of-mass whirl orbit
  17_Thermal_Settling.png — LPTN temperature convergence curve
"""

println("="^60)
println("  Julia ADORE — Transient Dynamics Plots")
println("="^60)
flush(stdout)

t_wall = time()

using ADORE
using Printf
using Statistics: mean
using LinearAlgebra: norm
using Plots;
ENV["GKSwstype"] = "100"
gr();

# --- Publication Quality Defaults ---
default(
    fontfamily="Computer Modern",
    linewidth=1.5,
    framestyle=:box,
    grid=true, gridalpha=0.3,
    titlefontsize=14, guidefontsize=12, tickfontsize=10, legendfontsize=10,
    dpi=300
)

# ──────────────────────────────────────────────────────────
# 1. Setup & Run Simulation
# ──────────────────────────────────────────────────────────
case_file = length(ARGS) >= 1 ? ARGS[1] : "Case4_transient.toml"
config_file = joinpath(@__DIR__, "inputs", case_file)
case_name = replace(case_file, ".toml" => "")

println("\n  Processing $case_name")
println("  Config: $config_file")

geom, mat, lub, trac, cage, config = load_simulation_config(config_file)

n_rpm = config.inner_race_speed * 30 / π
F_a = config.F_axial
Z = geom.n_balls

@printf("  Bearing: %d balls (pitch dia: %.1f mm)\n", Z, geom.d_m * 1000)
@printf("  Load: Fa=%.0f N, n=%.0f RPM\n", F_a, n_rpm)
@printf("  Duration: %.1f s (%.0f ms)\n", config.t_end, config.t_end * 1000)
flush(stdout)

println("\nRunning simulation (this may take several minutes)...")
flush(stdout)
result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)

elapsed = time() - t_wall
@printf("\n  Wall time: %.1f s\n", elapsed)
@printf("  ODE retcode: %s\n", result.retcode)
@printf("  Time points saved: %d\n", length(result.t))

# ──────────────────────────────────────────────────────────
# 2. Extract Data
# ──────────────────────────────────────────────────────────
t_dim = result.t
u_mat = result.u
scales = result.scales
n_t = length(t_dim)
p = result.params
to_rpm = 30.0 / π
t_ms = t_dim .* 1e3
t_s = t_dim

# ── Ball spin axis (body-frame ω) for Ball 1 ──
ball_idx = 1  # focus on Ball 1
bv_off = ADORE.ball_vel_offset(ball_idx, Z)
ω_x_all = u_mat[:, bv_off+3] .* scales.W   # axial body spin [rad/s]
ω_r_all = u_mat[:, bv_off+4] .* scales.W   # radial body spin
ω_θ_all = u_mat[:, bv_off+5] .* scales.W   # tangential body spin

# Spin axis pitch angle: angle between spin vector and axial direction
ω_mag = sqrt.(ω_x_all.^2 .+ ω_r_all.^2 .+ ω_θ_all.^2 .+ 1e-30)
β_spin = rad2deg.(atan.(sqrt.(ω_r_all.^2 .+ ω_θ_all.^2), ω_x_all))  # [°]

# ── Cage orbit ──
cage_poff = ADORE.cage_pos_offset(Z)
cage_y_um = u_mat[:, cage_poff+1] .* scales.L .* 1e6
cage_z_um = u_mat[:, cage_poff+2] .* scales.L .* 1e6
cage_ecc_um = sqrt.(cage_y_um.^2 .+ cage_z_um.^2)

cage_voff = ADORE.cage_vel_offset(Z)
cage_omega = u_mat[:, cage_voff+3] .* scales.W
cage_rpm = cage_omega .* to_rpm

# Theoretical cage speed
gamma_cos = (geom.d / geom.d_m) * cos(geom.alpha_0)
cage_theory_rpm = n_rpm * 0.5 * (1 - gamma_cos)
cage_err_pct = (cage_rpm .- cage_theory_rpm) ./ cage_theory_rpm .* 100

# ── Temperatures ──
th_off = ADORE.thermal_offset(Z)
T_ir_C = u_mat[:, th_off] .- 273.15
T_or_C = u_mat[:, th_off+1] .- 273.15
T_ball_C = u_mat[:, th_off+2] .- 273.15
T_oil_C = u_mat[:, th_off+3] .- 273.15

# ──────────────────────────────────────────────────────────
# 3. Generate Publication Plots
# ──────────────────────────────────────────────────────────
out_dir = joinpath(@__DIR__, "results", case_name)
mkpath(out_dir)
println("\n=== Generating Transient Plots → $out_dir ===")

# ── FIGURE 15: Spin Axis Nutation ──
# Left: ω_x vs ω_r phase portrait (colored by time)
# Right: Spin pitch angle β(t) time history

# Subsample for scatter coloring
n_sub = min(2000, n_t)
idx_sub = round.(Int, range(1, n_t, length=n_sub))

# Phase portrait
fig15a = scatter(ω_x_all[idx_sub] .* to_rpm, ω_r_all[idx_sub] .* to_rpm,
    marker_z=t_s[idx_sub], color=:viridis, ms=1.5, markerstrokewidth=0,
    colorbar_title="Time (s)", label="",
    xlabel="ω_x (RPM — axial)", ylabel="ω_r (RPM — radial)",
    title="Ball 1 Spin Axis Phase Portrait",
    aspect_ratio=:equal)
# Mark start and end
scatter!(fig15a, [ω_x_all[1]*to_rpm], [ω_r_all[1]*to_rpm],
    color=:lime, ms=6, markershape=:diamond, label="t=0 (QS init)")
scatter!(fig15a, [ω_x_all[end]*to_rpm], [ω_r_all[end]*to_rpm],
    color=:red, ms=6, markershape=:star5, label="t=end (steady)")

# Time history of pitch angle
fig15b = plot(t_s .* 1e3, β_spin,
    color=:darkblue, lw=1.2, label="",
    xlabel="Time (ms)", ylabel="Spin Axis Pitch β (°)",
    title="Spin Axis Nutation Settling")
# Mark the steady-state value
β_ss = mean(β_spin[max(1,n_t-100):n_t])
hline!(fig15b, [β_ss], color=:red, linestyle=:dash,
    label=@sprintf("Steady β = %.1f°", β_ss), lw=1.5)

# Inset zoom on first 200ms
t_zoom = 200.0  # ms
i_zoom = findfirst(t_ms .> t_zoom)
if i_zoom !== nothing && i_zoom > 20
    fig15c = plot(t_ms[1:i_zoom], β_spin[1:i_zoom],
        color=:darkblue, lw=1.5, label="",
        xlabel="Time (ms)", ylabel="β (°)",
        title="Early Nutation (0–$(Int(t_zoom)) ms)")
    hline!(fig15c, [β_ss], color=:red, linestyle=:dash, label="", lw=1)
    fig15 = plot(fig15a, plot(fig15b, fig15c, layout=(2,1)),
        layout=grid(1, 2, widths=[0.45, 0.55]), size=(1400, 600))
else
    fig15 = plot(fig15a, fig15b, layout=(1,2), size=(1400, 500))
end

savefig(fig15, joinpath(out_dir, "15_Spin_Nutation.png"))
println("  [OK] 15_Spin_Nutation.png")

# ── FIGURE 16: Cage Whirl Orbit ──
# Left: y_c vs z_c orbit colored by time
# Right top: eccentricity e_c(t), Right bottom: cage speed error

# Orbit plot
fig16a = scatter(cage_y_um[idx_sub], cage_z_um[idx_sub],
    marker_z=t_s[idx_sub], color=:plasma, ms=1.5, markerstrokewidth=0,
    colorbar_title="Time (s)", label="",
    xlabel="y_cage (μm)", ylabel="z_cage (μm)",
    title="Cage Center Whirl Orbit",
    aspect_ratio=:equal)
scatter!(fig16a, [cage_y_um[1]], [cage_z_um[1]],
    color=:lime, ms=6, markershape=:diamond, label="t=0")
scatter!(fig16a, [cage_y_um[end]], [cage_z_um[end]],
    color=:red, ms=6, markershape=:star5, label="t=end")

# Pilot clearance circle
θ_circ = range(0, 2π, length=100)
pilot_clr_um = cage.pilot_clearance * 1e6
plot!(fig16a, pilot_clr_um .* cos.(θ_circ), pilot_clr_um .* sin.(θ_circ),
    color=:red, linestyle=:dash, label="Pilot clearance", lw=1.5)

# Auto-scale orbit view
orbit_max = max(maximum(abs.(cage_y_um)), maximum(abs.(cage_z_um)), pilot_clr_um) * 1.3
plot!(fig16a, xlim=(-orbit_max, orbit_max), ylim=(-orbit_max, orbit_max))

# Eccentricity time history
fig16b = plot(t_ms, cage_ecc_um,
    color=:teal, lw=1.0, label="",
    xlabel="Time (ms)", ylabel="Cage Eccentricity (μm)",
    title="Cage Eccentricity vs Time")
hline!(fig16b, [pilot_clr_um], color=:red, linestyle=:dash,
    label="Pilot clearance", lw=1.5)

# Cage speed error
fig16c = plot(t_ms, cage_err_pct,
    color=:orange, lw=1.0, label="",
    xlabel="Time (ms)", ylabel="Cage Speed Error (%)",
    title="Cage Speed vs Harris Theory")
hline!(fig16c, [0.0], color=:black, linestyle=:dot, label="", lw=0.5)

fig16 = plot(fig16a, plot(fig16b, fig16c, layout=(2,1)),
    layout=grid(1, 2, widths=[0.45, 0.55]), size=(1400, 600))
savefig(fig16, joinpath(out_dir, "16_Cage_Whirl.png"))
println("  [OK] 16_Cage_Whirl.png")

# ── FIGURE 17: Thermal Settling ──
fig17 = plot(title="LPTN Temperature Response (100× accelerated thermal mass)",
    xlabel="Time (ms)", ylabel="Temperature (°C)",
    legend=:right, size=(1000, 500))
plot!(fig17, t_ms, T_ir_C, color=:steelblue, lw=2, label="Inner Race")
plot!(fig17, t_ms, T_ball_C, color=:darkred, lw=2, label="Ball (lumped)")
plot!(fig17, t_ms, T_or_C, color=:orange, lw=2, label="Outer Race")
plot!(fig17, t_ms, T_oil_C, color=:darkgreen, lw=2, linestyle=:dash, label="Lubricant")

# Annotate final temperatures
@printf("  Final T: IR=%.1f°C, OR=%.1f°C, Ball=%.1f°C, Oil=%.1f°C\n",
    T_ir_C[end], T_or_C[end], T_ball_C[end], T_oil_C[end])

savefig(fig17, joinpath(out_dir, "17_Thermal_Settling.png"))
println("  [OK] 17_Thermal_Settling.png")

# ── Summary ──
println("\n" * "="^60)
@printf("  Ball 1 spin axis pitch (steady): β = %.1f°\n", β_ss)
@printf("  Cage eccentricity (final): %.1f μm (limit: %.1f μm)\n",
    cage_ecc_um[end], pilot_clr_um)
@printf("  Cage speed error (final): %.3f%%\n", cage_err_pct[end])
@printf("  Final cage speed: %.1f RPM (theory: %.1f)\n", cage_rpm[end], cage_theory_rpm)
@printf("  Wall time: %.1f s\n", time() - t_wall)
println("="^60)
println("  TRANSIENT PLOTS COMPLETE → $out_dir")
println("="^60)
