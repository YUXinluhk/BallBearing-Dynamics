#!/usr/bin/env julia
"""
run_postprocess.jl — Run simulation and generate postprocessing plots.

Matches output format of ADORE_V2/run_dynamic_postprocess.py.
Plots: 01_Race_Displacement, 02_Cage_Dynamics, 03_Contact_Loads, 04_Ball_Kinematics
       07_All_Ball_Loads, 08_All_Ball_Forces, 09_All_Ball_Moments
       10_Heat_Generation, 11_All_Ball_Heat, 12_Heat_Pie, 13_Heat_Per_Ball
"""

println("="^60)
println("  Julia ADORE — Dynamic Postprocessing")
println("="^60)
flush(stdout)

t_wall = time()

using ADORE
using Printf
using Statistics: mean
using LinearAlgebra: norm
using Plots;
ENV["GKSwstype"] = "100"
gr();  # GR backend

# ──────────────────────────────────────────────────────────
# 1. Setup & Run Simulation
# ──────────────────────────────────────────────────────────
geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

n_rpm = 10000.0
F_a = 2000.0
omega = n_rpm * π / 30

config = SimulationConfig(
    t_end=0.005,
    dt_output=50e-6,
    inner_race_speed=omega,
    F_axial=F_a,
    t_ramp_end=0.005,
    zeta=0.10,
    c_structural=0.0,
    mu_spin=0.06,
    integrator=IntegratorConfig(
        rtol=1e-4, atol=1e-7, h_max=1e-4, max_steps=10_000_000,
    ),
)

@printf("  Bearing: 7210B (%d balls)\n", geom.n_balls)
@printf("  Load: Fa=%.0f N, n=%.0f RPM\n", F_a, n_rpm)
@printf("  Duration: %.0f ms\n", config.t_end * 1000)
flush(stdout)

println("\nRunning simulation...")
flush(stdout)
result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)

elapsed = time() - t_wall
@printf("\n  Wall time: %.1f s\n", elapsed)
@printf("  ODE retcode: %s\n", result.retcode)
@printf("  Time points: %d\n", length(result.t))

# ──────────────────────────────────────────────────────────
# 2. Extract Kinematics
# ──────────────────────────────────────────────────────────
t_dim = result.t     # dimensional time [s]
u_mat = result.u     # (n_t, n_state) nondim
Z = geom.n_balls
scales = result.scales
n_t = length(t_dim)
p = result.params     # flat parameter vector from kernel

t_ms = t_dim .* 1e3
to_rpm = 30.0 / π

# ── Inner race (dimensional) ──
x_ir_um = u_mat[:, 1] .* scales.L .* 1e6   # μm
y_ir_um = u_mat[:, 2] .* scales.L .* 1e6
z_ir_um = u_mat[:, 3] .* scales.L .* 1e6

ir_v_off = ADORE.ir_vel_offset(Z)
vx_ir = u_mat[:, ir_v_off] .* scales.V

# ── Cage ──
cage_voff = ADORE.cage_vel_offset(Z)
cage_omega = u_mat[:, cage_voff+3] .* scales.W   # rad/s
cage_rpm = cage_omega .* to_rpm

cage_poff = ADORE.cage_pos_offset(Z)
cage_y_um = u_mat[:, cage_poff+1] .* scales.L .* 1e6
cage_z_um = u_mat[:, cage_poff+2] .* scales.L .* 1e6

# Theoretical cage speed
gamma_cos = (geom.d / geom.d_m) * cos(geom.alpha_0)
cage_theory_rpm = n_rpm * 0.5 * (1 - gamma_cos)

# ── Per-ball kinematics ──
orbit_rpm = zeros(n_t, Z)
spin_total_rpm = zeros(n_t, Z)
spin_roll_i_rpm = zeros(n_t, Z)
spin_roll_o_rpm = zeros(n_t, Z)
spin_norm_i_rpm = zeros(n_t, Z)
spin_norm_o_rpm = zeros(n_t, Z)
spin_y_rpm = zeros(n_t, Z)   # side-slip (ω_y body)

for j in 1:Z
    bv_off = ADORE.ball_vel_offset(j, Z)
    orbit_rpm[:, j] .= u_mat[:, bv_off+2] .* scales.W .* to_rpm  # θ̇ → RPM

    for i in 1:n_t
        ωx = u_mat[i, bv_off+3] * scales.W
        ωy = u_mat[i, bv_off+4] * scales.W
        ωz = u_mat[i, bv_off+5] * scales.W

        # Absolute Total Rotational Speed in X-Z cross-section
        spin_total_rpm[i, j] = sqrt(ωx^2 + ωz^2) * to_rpm

        # Fetch dynamic contact angles computed later (or we can use nominal for now, 
        # but to be perfectly accurate we will compute them inline first)
        bp = ADORE.ball_pos_offset(j)
        x_b, r_b, θ_b = u_mat[i, bp], u_mat[i, bp+1], u_mat[i, bp+2]
        sθ, cθ = sincos(θ_b)
        x_ir_s, y_ir_s, z_ir_s = u_mat[i, 1], u_mat[i, 2], u_mat[i, 3]

        # Inner GC
        r_gc_i = p[ADORE.P_DI] / 2 + y_ir_s * (-sθ) + z_ir_s * cθ
        x_gc_i = p[ADORE.P_X_GI0] + x_ir_s
        α_i = atan(x_gc_i - x_b, r_gc_i - r_b)

        # Outer GC
        r_gc_o = p[ADORE.P_DO] / 2
        x_gc_o = p[ADORE.P_X_GO0]
        α_o = atan(-(x_gc_o - x_b), -(r_gc_o - r_b))

        # True rolling vectors
        ω_roll_i = ωx * cos(α_i) - ωz * sin(α_i)
        ω_roll_o = ωx * cos(α_o) - ωz * sin(α_o)

        # True spin vectors (normal to contact)
        ω_spin_i = ωx * sin(α_i) + ωz * cos(α_i)
        ω_spin_o = ωx * sin(α_o) + ωz * cos(α_o)

        spin_roll_i_rpm[i, j] = ω_roll_i * to_rpm
        spin_roll_o_rpm[i, j] = ω_roll_o * to_rpm
        spin_norm_i_rpm[i, j] = ω_spin_i * to_rpm
        spin_norm_o_rpm[i, j] = ω_spin_o * to_rpm
        spin_y_rpm[i, j] = ωy * to_rpm             # Lateral pitch
    end
end
mean_orbit_rpm = vec(mean(orbit_rpm, dims=2))

# ── Contact geometry (from kernel params) ──
D_i_star = p[ADORE.P_DI]
D_o_star = p[ADORE.P_DO]
drb_i = p[ADORE.P_DRB_I]
drb_o = p[ADORE.P_DRB_O]
x_gi0_star = p[ADORE.P_X_GI0]
x_go0_star = p[ADORE.P_X_GO0]
Y_i = p[ADORE.P_YI]
Y_o = p[ADORE.P_YO]
Q_scale = scales.Q

Q_inner = zeros(n_t, Z)
Q_outer = zeros(n_t, Z)
alpha_inner = zeros(n_t, Z)
alpha_outer = zeros(n_t, Z)

for i in 1:n_t
    x_ir_s = u_mat[i, 1]
    y_ir_s = u_mat[i, 2]
    z_ir_s = u_mat[i, 3]

    for j in 1:Z
        bp = ADORE.ball_pos_offset(j)
        x_b = u_mat[i, bp]
        r_b = u_mat[i, bp+1]
        θ_b = u_mat[i, bp+2]

        sθ, cθ = sincos(θ_b)

        # Groove center inner (same as kernel)
        r_gc_i = D_i_star / 2 + y_ir_s * (-sθ) + z_ir_s * cθ
        x_gc_i = x_gi0_star + x_ir_s

        # Groove center outer (fixed)
        r_gc_o = D_o_star / 2
        x_gc_o = x_go0_star

        # Vectors from groove center to ball center
        dx_i = x_gc_i - x_b
        dr_i = r_gc_i - r_b
        dx_o = x_gc_o - x_b
        dr_o = r_gc_o - r_b

        # Distance
        L_i = sqrt(dx_i^2 + dr_i^2 + 1e-30)
        L_o = sqrt(dx_o^2 + dr_o^2 + 1e-30)

        # Penetration
        δ_i = max(0.0, L_i - drb_i)
        δ_o = max(0.0, L_o - drb_o)

        # Hertz loads (nondim → dimensional N)
        Q_inner[i, j] = Y_i * δ_i^1.5 * Q_scale
        Q_outer[i, j] = Y_o * δ_o^1.5 * Q_scale

        # Contact angles (degrees)
        alpha_inner[i, j] = atand(dx_i, dr_i)
        alpha_outer[i, j] = atand(-dx_o, -dr_o)
    end
end

# Ball azimuthal positions (degrees)
theta_deg = zeros(n_t, Z)
for j in 1:Z
    bp = ADORE.ball_pos_offset(j)
    theta_deg[:, j] .= mod.(rad2deg.(u_mat[:, bp+2]), 360.0)
end

# ──────────────────────────────────────────────────────────
# 3. Generate Plots
# ──────────────────────────────────────────────────────────
out_dir = joinpath(@__DIR__, "results")
mkpath(out_dir)

println("\n=== Generating Plots → $out_dir ===")

# ── Plot 01: Race Displacement ──
p1 = plot(t_ms, x_ir_um, lw=1.5, color=:black, legend=false,
    title="Inner Race Axial Displacement", ylabel="Axial X (μm)")
p2 = plot(t_ms, y_ir_um, label="Radial Y", color=:steelblue, lw=1.2, alpha=0.8)
plot!(p2, t_ms, z_ir_um, label="Radial Z", color=:red, lw=1.2, alpha=0.8,
    title="Inner Race Radial Displacement", xlabel="Time (ms)", ylabel="Radial (μm)")
fig01 = plot(p1, p2, layout=(2, 1), size=(800, 600))
savefig(fig01, joinpath(out_dir, "01_Race_Displacement.png"))
println("  [OK] 01_Race_Displacement.png")

# ── Plot 02: Cage Dynamics ──
p3 = plot(t_ms, cage_rpm, label="Actual Cage Speed", color=:orange, lw=2)
hline!(p3, [cage_theory_rpm], color=:black, linestyle=:dash, label="Theoretical",
    title="Cage Speed vs Theoretical", xlabel="Time (ms)", ylabel="Speed (RPM)")

p4 = plot(cage_y_um, cage_z_um, color=:teal, alpha=0.7, lw=0.8, label="",
    aspect_ratio=:equal, title="Cage Whirl Orbit", xlabel="Y (μm)", ylabel="Z (μm)")
scatter!(p4, [cage_y_um[end]], [cage_z_um[end]], color=:red, label="End", ms=4)
# Draw pilot clearance circle
θ_circ = range(0, 2π, length=100)
pilot_um = 50.0
plot!(p4, pilot_um .* cos.(θ_circ), pilot_um .* sin.(θ_circ),
    color=:red, linestyle=:dash, label="Pilot Limit", xlim=(-80, 80), ylim=(-80, 80))

fig02 = plot(p3, p4, layout=(1, 2), size=(1200, 450))
savefig(fig02, joinpath(out_dir, "02_Cage_Dynamics.png"))
println("  [OK] 02_Cage_Dynamics.png")

# ── Plot 03: Contact Loads ──
# Find heaviest loaded ball
h_idx = argmax(vec(mean(Q_outer[max(1, n_t - 50):n_t, :], dims=1)))

p5 = plot(t_ms, Q_inner[:, h_idx], label="Inner (Ball $h_idx)", color=:steelblue, lw=1.5)
plot!(p5, t_ms, Q_outer[:, h_idx], label="Outer (Ball $h_idx)", color=:red, lw=1.5, linestyle=:dash,
    title="Transient Contact Loads (Heaviest Ball)", xlabel="Time (ms)", ylabel="Contact Load (N)")

p6 = bar(0:Z-1, [Q_inner[end, :] Q_outer[end, :]], label=["Inner" "Outer"],
    color=[:steelblue :red], alpha=0.8,
    title="Load Distribution (Final Step)", xlabel="Ball Number", ylabel="Load (N)")

# Polar load zone
θ_fin = deg2rad.(theta_deg[end, :])
Q_fin = Q_outer[end, :]
idx = sortperm(θ_fin)
θ_sorted = [θ_fin[idx]; θ_fin[idx[1]]]
Q_sorted = [Q_fin[idx]; Q_fin[idx[1]]]
p7 = plot(θ_sorted, Q_sorted, proj=:polar, color=:darkred, marker=:circle, ms=3,
    fill=true, fillalpha=0.2, fillcolor=:red,
    title="Outer Race Load Zone (Polar)")

fig03 = plot(p5, plot(p6, p7, layout=(1, 2)), layout=@layout([a; b]),
    size=(1200, 800))
savefig(fig03, joinpath(out_dir, "03_Contact_Loads.png"))
println("  [OK] 03_Contact_Loads.png")

# ── Plot 04: Ball Kinematics ──
pp1 = plot(t_ms, spin_total_rpm[:, 1], color=:black, lw=1.5,
    title="Total Ball Rotational Speed (ωx, ωz magnitude)", ylabel="Total (RPM)", label="Ball 1 Total")
plot!(pp1, t_ms, mean_orbit_rpm, color=:gray, lw=1.0, linestyle=:dash, label="Mean Orbital Speed")

pp2 = plot(t_ms, spin_roll_i_rpm[:, 1], color=:green, lw=1.2,
    title="True Rolling Speed (relative to normal)", ylabel="Roll (RPM)", label="Ball 1 Inner Roll")
plot!(pp2, t_ms, spin_roll_o_rpm[:, 1], color=:darkgreen, lw=1.2, linestyle=:dash, label="Ball 1 Outer Roll")

pp3 = plot(t_ms, spin_norm_i_rpm[:, 1], color=:purple, lw=1.2,
    title="True Gyroscopic Spin (along normal)", ylabel="Spin (RPM)", label="Ball 1 Inner Spin")
plot!(pp3, t_ms, spin_norm_o_rpm[:, 1], color=:indigo, lw=1.2, linestyle=:dash, label="Ball 1 Outer Spin")

pp4 = plot(t_ms, spin_y_rpm[:, 1], color=:red, lw=1.2,
    title="Transverse Pitch / Side-Slip (ωy axis)", ylabel="Pitch (RPM)", label="Ball 1 Pitch Slip")

pp5 = plot(t_ms, alpha_inner[:, 1], label="Inner Contact Angle", color=:steelblue, lw=1.2)
plot!(pp5, t_ms, alpha_outer[:, 1], label="Outer Contact Angle", color=:red, lw=1.2, linestyle=:dash)
hline!(pp5, [rad2deg(geom.alpha_0)], color=:black, linestyle=:dot, label="Nominal Free Angle",
    title="Dynamic Contact Angles (Ball 1)", xlabel="Time (ms)", ylabel="Angle (°)")

fig04 = plot(pp1, pp2, pp3, pp4, pp5, layout=(5, 1), size=(800, 1500))
savefig(fig04, joinpath(out_dir, "04_Ball_Kinematics.png"))
println("  [OK] 04_Ball_Kinematics.png")

# ── Summary ──
@printf("\n  Theory cage speed: %.1f RPM\n", cage_theory_rpm)
@printf("  Final cage speed:  %.1f RPM (%.1f%%)\n", cage_rpm[end], cage_rpm[end] / cage_theory_rpm * 100)
@printf("  Final Q_inner(Ball 1): %.1f N\n", Q_inner[end, 1])
@printf("  Final α_inner(Ball 1): %.1f°, α_outer(Ball 1): %.1f°\n", alpha_inner[end, 1], alpha_outer[end, 1])

# ══════════════════════════════════════════════════════════
# PLOTS 07–13: All-Ball Forces, Moments & Heat Generation
# ══════════════════════════════════════════════════════════

println("\n=== Computing Field Outputs ===")
flush(stdout)

fo_matrix = compute_field_outputs(result)   # (n_t, Z × N_FIELD_PER_BALL)

# Reshape to per-ball arrays: fo_ball[i, j, field]
fo_ball = zeros(n_t, Z, ADORE.N_FIELD_PER_BALL)
for j in 1:Z
    base = (j - 1) * ADORE.N_FIELD_PER_BALL
    for f in 1:ADORE.N_FIELD_PER_BALL
        fo_ball[:, j, f] .= fo_matrix[:, base+f]
    end
end

# Field outputs are now in DIMENSIONAL units (fixed in field_output_kernel)
# Extract directly — no scaling needed

# Contact loads [N] and angles [rad]
fo_Qi = fo_ball[:, :, ADORE.FO_Q_I]
fo_Qo = fo_ball[:, :, ADORE.FO_Q_O]
fo_alpha_i = fo_ball[:, :, ADORE.FO_ALPHA_I]
fo_alpha_o = fo_ball[:, :, ADORE.FO_ALPHA_O]

# Forces [N]
fo_F_trac_i = fo_ball[:, :, ADORE.FO_F_TRAC_I]
fo_F_trac_o = fo_ball[:, :, ADORE.FO_F_TRAC_O]
fo_F_drag = fo_ball[:, :, ADORE.FO_F_DRAG]
fo_F_pocket = fo_ball[:, :, ADORE.FO_F_POCKET]

# Moments [N·m]
fo_M_spin_i = fo_ball[:, :, ADORE.FO_M_SPIN_I]
fo_M_spin_o = fo_ball[:, :, ADORE.FO_M_SPIN_O]

# Heat sources [W]
fo_H_slide_i = fo_ball[:, :, ADORE.FO_H_SLIDE_I]
fo_H_slide_o = fo_ball[:, :, ADORE.FO_H_SLIDE_O]
fo_H_spin_i = fo_ball[:, :, ADORE.FO_H_SPIN_I]
fo_H_spin_o = fo_ball[:, :, ADORE.FO_H_SPIN_O]
fo_H_drag = fo_ball[:, :, ADORE.FO_H_DRAG]
fo_H_churn = fo_ball[:, :, ADORE.FO_H_CHURN]
fo_H_total = fo_H_slide_i .+ fo_H_slide_o .+ fo_H_spin_i .+ fo_H_spin_o .+ fo_H_drag .+ fo_H_churn

# Angular velocities [rad/s]
fo_omega_x = fo_ball[:, :, ADORE.FO_OMEGA_X]
fo_omega_y = fo_ball[:, :, ADORE.FO_OMEGA_Y]
fo_omega_z = fo_ball[:, :, ADORE.FO_OMEGA_Z]
fo_w_spin_i = fo_ball[:, :, ADORE.FO_W_SPIN_I]
fo_w_spin_o = fo_ball[:, :, ADORE.FO_W_SPIN_O]
fo_theta_dot = fo_ball[:, :, ADORE.FO_THETA_DOT]

# Gyroscopic moment: M_g = I_ball * ω_spin * ω_precession [N·m]
I_ball = ball_inertia(geom)  # kg·m²
fo_omega_mag = sqrt.(fo_omega_x .^ 2 .+ fo_omega_y .^ 2 .+ fo_omega_z .^ 2)
fo_M_gyro = I_ball .* fo_omega_mag .* abs.(fo_theta_dot)
D_ball = geom.d   # ball diameter [m]
fo_MgD = fo_M_gyro ./ D_ball   # M_g/D [N]

@printf("  Field outputs: %d time steps × %d balls\n", n_t, Z)
flush(stdout)

# Define color palette for balls (tab20 equivalent, hex strings)
ball_colors = [
    "#1f77b4", "#aec7e8",   # B0, B1
    "#ff7f0e", "#ffbb78",   # B2, B3
    "#2ca02c", "#98df8a",   # B4, B5
    "#d62728", "#ff9896",   # B6, B7
    "#9467bd", "#c5b0d5",   # B8, B9
    "#8c564b", "#c49c94",   # B10, B11
    "#e377c2", "#f7b6d2",   # B12, B13
    "#7f7f7f", "#c7c7c7",   # B14, B15 (spare)
    "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5",
]

println("\n=== Generating Plots 07–13 ===")
flush(stdout)

# ── Plot 07: All Ball Loads ──
fig07_p1 = plot(title="Qᵢ — Inner Race Contact Load (N)", ylabel="", grid=true, gridalpha=0.3)
fig07_p2 = plot(title="Qₒ — Outer Race Contact Load (N)", ylabel="", grid=true, gridalpha=0.3)
fig07_p3 = plot(title="αᵢ — Inner Contact Angle (°)", ylabel="", xlabel="Time (ms)", grid=true, gridalpha=0.3)
fig07_p4 = plot(title="αₒ — Outer Contact Angle (°)", ylabel="", xlabel="Time (ms)", grid=true, gridalpha=0.3)

for j in 1:Z
    c = ball_colors[min(j, length(ball_colors))]
    lbl = "B$(j-1)"
    plot!(fig07_p1, t_ms, fo_Qi[:, j], color=c, label=lbl, linewidth=0.5, legend=false)
    plot!(fig07_p2, t_ms, fo_Qo[:, j], color=c, label=lbl, linewidth=0.5, legend=false)
    plot!(fig07_p3, t_ms, rad2deg.(fo_alpha_i[:, j]), color=c, label=lbl, linewidth=0.5, legend=false)
    plot!(fig07_p4, t_ms, rad2deg.(fo_alpha_o[:, j]), color=c, label=lbl, linewidth=0.5, legend=false)
end

fig07 = plot(fig07_p1, fig07_p2, fig07_p3, fig07_p4, layout=(2, 2), size=(1600, 900),
    plot_title="All Balls — Contact Load History")
savefig(fig07, joinpath(out_dir, "07_All_Ball_Loads.png"))
println("  [OK] 07_All_Ball_Loads.png")

# ── Plot 08: All Ball Forces ──
# Normal force sum (axial/radial) — using contact angle projections
fig08_p1 = plot(title="Normal Force Sum — Axial (N)", ylabel="", grid=true, gridalpha=0.3)
fig08_p2 = plot(title="Normal Force Sum — Radial (N)", ylabel="", grid=true, gridalpha=0.3)
fig08_p3 = plot(title="Traction Force — Axial (N)", ylabel="", grid=true, gridalpha=0.3)
fig08_p4 = plot(title="Traction Force — Tangential (N)", ylabel="", grid=true, gridalpha=0.3)
fig08_p5 = plot(title="Cage Pocket Force — Axial (N)", ylabel="", xlabel="Time (ms)", grid=true, gridalpha=0.3)
fig08_p6 = plot(title="Churning Drag — Tangential (N)", ylabel="", xlabel="Time (ms)", grid=true, gridalpha=0.3)

for j in 1:Z
    c = ball_colors[min(j, length(ball_colors))]
    lbl = "B$(j-1)"

    # Normal force: Q_i * sin(α_i) + Q_o * sin(α_o) for axial
    F_ax = fo_Qi[:, j] .* sin.(fo_alpha_i[:, j]) .- fo_Qo[:, j] .* sin.(fo_alpha_o[:, j])
    F_rad = fo_Qi[:, j] .* cos.(fo_alpha_i[:, j]) .- fo_Qo[:, j] .* cos.(fo_alpha_o[:, j])
    plot!(fig08_p1, t_ms, F_ax, color=c, linewidth=0.5, label=lbl, legend=false)
    plot!(fig08_p2, t_ms, F_rad, color=c, linewidth=0.5, label=lbl, legend=false)

    # Traction force (approx: use F_trac magnitude for axial/tangential approximation)
    plot!(fig08_p3, t_ms, fo_F_trac_i[:, j], color=c, linewidth=0.5, label=lbl, legend=false)
    plot!(fig08_p4, t_ms, fo_F_trac_o[:, j], color=c, linewidth=0.5, label=lbl, legend=false)

    # Cage pocket & drag
    plot!(fig08_p5, t_ms, fo_F_pocket[:, j], color=c, linewidth=0.5, label=lbl, legend=false)
    plot!(fig08_p6, t_ms, fo_F_drag[:, j], color=c, linewidth=0.5, label=lbl, legend=false)
end

fig08 = plot(fig08_p1, fig08_p2, fig08_p3, fig08_p4, fig08_p5, fig08_p6,
    layout=(3, 2), size=(1600, 1200), plot_title="All Balls — Force Components")
savefig(fig08, joinpath(out_dir, "08_All_Ball_Forces.png"))
println("  [OK] 08_All_Ball_Forces.png")

# ── Plot 09: All Ball Moments ──
fig09_p1 = plot(title="M_spin,inner (N·mm)", ylabel="", grid=true, gridalpha=0.3)
fig09_p2 = plot(title="M_spin,outer (N·mm)", ylabel="", grid=true, gridalpha=0.3)
fig09_p3 = plot(title="M_gyro (N·mm)", ylabel="", xlabel="Time (ms)", grid=true, gridalpha=0.3)
fig09_p4 = plot(title="M_g/D (N)", ylabel="", xlabel="Time (ms)", grid=true, gridalpha=0.3)

for j in 1:Z
    c = ball_colors[min(j, length(ball_colors))]
    lbl = "B$(j-1)"
    plot!(fig09_p1, t_ms, fo_M_spin_i[:, j] .* 1e3, color=c, linewidth=0.5, label=lbl, legend=false)
    plot!(fig09_p2, t_ms, fo_M_spin_o[:, j] .* 1e3, color=c, linewidth=0.5, label=lbl, legend=false)
    plot!(fig09_p3, t_ms, fo_M_gyro[:, j] .* 1e3, color=c, linewidth=0.5, label=lbl, legend=false)
    plot!(fig09_p4, t_ms, fo_MgD[:, j], color=c, linewidth=0.5, label=lbl, legend=false)
end

fig09 = plot(fig09_p1, fig09_p2, fig09_p3, fig09_p4, layout=(2, 2), size=(1600, 900),
    plot_title="All Balls — Moment Components")
savefig(fig09, joinpath(out_dir, "09_All_Ball_Moments.png"))
println("  [OK] 09_All_Ball_Moments.png")

# ── Plot 10: Ball 0 — Heat Generation Breakdown ──
b0 = 1  # Julia 1-indexed, Ball 0

heat_sources = [
    (fo_H_slide_i[:, b0], "Inner Sliding", "#1f77b4"),
    (fo_H_slide_o[:, b0], "Outer Sliding", "#d62728"),
    (fo_H_spin_i[:, b0], "Inner Spin", "#2ca02c"),
    (fo_H_spin_o[:, b0], "Outer Spin", "#ff7f0e"),
    (fo_H_drag[:, b0], "Churning Drag", "#9467bd"),
    (fo_H_churn[:, b0], "Churning Rot.", "#8c564b"),
]

# Top panel: Stacked area
fig10_top = plot(title="Heat Sources — Stacked Area (W)", ylabel="Power (W)", grid=true, gridalpha=0.3)
y_stack = zeros(n_t)
for (vals, label, color) in heat_sources
    global y_stack
    y_new = y_stack .+ vals
    plot!(fig10_top, t_ms, y_new, fillrange=y_stack, fillalpha=0.7, fillcolor=color,
        color=color, label=label, linewidth=0.5)
    y_stack = copy(y_new)
end

# Bottom panel: Individual lines
fig10_bot = plot(title="Heat Sources — Individual (W)", xlabel="Time (ms)",
    ylabel="Power (W)", grid=true, gridalpha=0.3)
for (vals, label, color) in heat_sources
    plot!(fig10_bot, t_ms, vals, label=label, color=color, linewidth=1.2)
end
plot!(fig10_bot, t_ms, fo_H_total[:, b0], label="Total", color=:black, linewidth=2, alpha=0.8)

fig10 = plot(fig10_top, fig10_bot, layout=(2, 1), size=(1200, 800),
    plot_title="Ball 0 — Friction Heat Generation Breakdown")
savefig(fig10, joinpath(out_dir, "10_Heat_Generation.png"))
println("  [OK] 10_Heat_Generation.png")

# ── Plot 11: All Ball Total Heat ──
fig11_top = plot(title="Per-Ball Total Heat (W)", ylabel="Power (W)", grid=true, gridalpha=0.3)
for j in 1:Z
    c = ball_colors[min(j, length(ball_colors))]
    plot!(fig11_top, t_ms, fo_H_total[:, j], color=c, label="B$(j-1)", linewidth=0.7, legend=false)
end

H_bearing_total = vec(sum(fo_H_total, dims=2))
fig11_bot = plot(t_ms, H_bearing_total, color=:black, linewidth=2, legend=false,
    title="Bearing Total Heat Generation (W)", xlabel="Time (ms)", ylabel="Power (W)",
    grid=true, gridalpha=0.3)

fig11 = plot(fig11_top, fig11_bot, layout=(2, 1), size=(1200, 700),
    plot_title="All Balls — Total Heat Generation")
savefig(fig11, joinpath(out_dir, "11_All_Ball_Heat.png"))
println("  [OK] 11_All_Ball_Heat.png")

# ── Plot 12: Steady-State Heat Pie Chart ──
n_steady = max(1, n_t ÷ 5)   # last 20%
ss_range = max(1, n_t - n_steady + 1):n_t

source_labels_data = [
    ("Inner Sliding", fo_H_slide_i, "#1f77b4"),
    ("Outer Sliding", fo_H_slide_o, "#d62728"),
    ("Inner Spin", fo_H_spin_i, "#2ca02c"),
    ("Outer Spin", fo_H_spin_o, "#ff7f0e"),
    ("Churning Drag", fo_H_drag, "#9467bd"),
    ("Churning Rot.", fo_H_churn, "#8c564b"),
]

ss_means = Float64[]
ss_labels = String[]
ss_colors = String[]
for (label, data, color) in source_labels_data
    val = mean(sum(data[ss_range, :], dims=2))
    push!(ss_means, val)
    push!(ss_labels, label)
    push!(ss_colors, color)
end

total_heat = sum(ss_means)

if total_heat > 1e-12
    # Filter out negligible sources (< 0.1%)
    keep = ss_means .> 0.001 * total_heat
    pie_vals = ss_means[keep]
    pie_labels_f = ss_labels[keep]
    pie_colors_f = ss_colors[keep]
    pie_pcts = pie_vals ./ sum(pie_vals) .* 100

    # Manual pie chart using Shape arcs (Plots.jl has no built-in pie)
    fig12 = plot(aspect_ratio=:equal, legend=:outertopright, grid=false,
        axis=false, ticks=false, size=(800, 700),
        title=@sprintf("Steady-State Heat Source Breakdown\n(Total: %.2f W, last %d steps)",
            total_heat, n_steady))

    start_angle = π / 2
    for k in eachindex(pie_vals)
        global start_angle
        frac = pie_vals[k] / sum(pie_vals)
        end_angle = start_angle - frac * 2π
        n_arc = max(3, Int(ceil(frac * 100)))
        θ_range = range(start_angle, end_angle, length=n_arc)
        xs = [0.0; cos.(θ_range)]
        ys = [0.0; sin.(θ_range)]
        plot!(fig12, Shape(xs, ys), fillcolor=pie_colors_f[k], fillalpha=0.85,
            linecolor=:white, linewidth=1.5,
            label=@sprintf("%s (%.1f%%)", pie_labels_f[k], pie_pcts[k]))
        # Label at mid-angle
        mid_θ = (start_angle + end_angle) / 2
        if frac > 0.05
            annotate!(fig12, 0.65 * cos(mid_θ), 0.65 * sin(mid_θ),
                text(@sprintf("%.1f%%", pie_pcts[k]), :center, 8, :white))
        end
        start_angle = end_angle
    end

    savefig(fig12, joinpath(out_dir, "12_Heat_Pie.png"))
    println("  [OK] 12_Heat_Pie.png")
else
    println("  [SKIP] 12_Heat_Pie.png — total heat ≈ 0")
end

# ── Plot 13: Steady-State Per-Ball Heat Bar Chart ──
balls_x = collect(0:Z-1)
bar_bottom = zeros(Z)

fig13 = plot(title="Per-Ball Heat Distribution (Steady State)", xlabel="Ball Number",
    ylabel="Heat Generation (W)", grid=true, gridalpha=0.3, legend=:topright, size=(1200, 500))

for (label, data, color) in source_labels_data
    vals = vec(mean(data[ss_range, :], dims=1))
    bar!(fig13, balls_x, vals, bar_position=:stack, label=label, color=color, alpha=0.85)
end

fig13 = plot(fig13, plot_title="Steady-State Heat by Ball and Source")
savefig(fig13, joinpath(out_dir, "13_Heat_Per_Ball.png"))
println("  [OK] 13_Heat_Per_Ball.png")

# ── Heat Summary ──
println("\n  === Friction Heat Summary (Steady-State, last $n_steady steps) ===")
for (label, val) in zip(ss_labels, ss_means)
    @printf("    %-20s: %8.3f W\n", label, val)
end
@printf("    %-20s: %8.3f W\n", "TOTAL", total_heat)
@printf("    Per-ball mean     : %8.3f W\n", total_heat / Z)
ss_per_ball = vec(mean(fo_H_total[ss_range, :], dims=1))
@printf("    Per-ball max      : %8.3f W\n", maximum(ss_per_ball))
@printf("    Per-ball min      : %8.3f W\n", minimum(ss_per_ball))

println("\n" * "="^60)
println("  POSTPROCESSING COMPLETE → $out_dir")
println("="^60)

