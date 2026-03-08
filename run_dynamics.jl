#!/usr/bin/env julia
"""
ADORE.jl — 7210B 1-Second Transient Simulation
Solver: FBDF implicit (BDF)
"""

println("="^60)
println("  ADORE.jl — 1s Physics Simulation (FBDF)")
println("="^60)
flush(stdout)

t_wall_start = time()

using ADORE
using Dates
using Printf
using Statistics: mean, std
using LinearAlgebra: norm

# ── Bearing setup ──
geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

n_rpm = 10000.0
F_a = 2000.0
t_end = 0.05

config = SimulationConfig(
    t_end=t_end,
    dt_output=1e-3,
    inner_race_speed=n_rpm * π / 30,
    F_axial=F_a,
    t_ramp_end=5e-3,
    zeta=0.03,
    c_structural=20.0,
    mu_spin=0.06,
    integrator=IntegratorConfig(
        rtol=1e-5,
        atol=1e-7,
        h_max=1e-4,
        max_steps=50_000_000,
    ),
)

ω_ir = n_rpm * π / 30
@printf("Bearing:    7210B (%d balls)\n", geom.n_balls)
@printf("Load:       F_a=%.0f N, n=%.0f RPM\n", F_a, n_rpm)
@printf("Duration:   %.3f s\n", t_end)
@printf("Solver:     FBDF + FiniteDiff Jacobian\n")
flush(stdout)

# ── Run simulation ──
println("\nStarting at ", Dates.format(Dates.now(), "HH:MM:SS"), "...")
flush(stdout)

result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)

t_wall = time() - t_wall_start

# ── Report ──
println("\n" * "="^60)
println("  RESULTS (", Dates.format(Dates.now(), "HH:MM:SS"), ")")
println("="^60)
@printf("  ODE retcode:  %s\n", result.retcode)
@printf("  Time points:  %d\n", length(result.t))
@printf("  Wall clock:   %.1f s (%.2f min)\n", t_wall, t_wall / 60)

t_dim = result.t
u_mat = result.u   # (n_time × n_state) matrix
Z = geom.n_balls
scales = result.scales
n_t = length(t_dim)

@printf("  u_mat size:   %s\n", string(size(u_mat)))
has_nan = any(isnan.(u_mat))
@printf("  NaN in sol:   %s\n\n", has_nan)

if !has_nan && n_t > 10
    try
        # ── Inner Race Axial Displacement ──
        println("--- Inner Race ---")
        # u_mat[time_idx, state_idx] — state idx 1 = IR axial position
        x_ir_0 = u_mat[1, 1] * scales.L
        x_ir_end = u_mat[end, 1] * scales.L
        x_ir_min = minimum(u_mat[:, 1]) * scales.L
        x_ir_max = maximum(u_mat[:, 1]) * scales.L
        @printf("  x_IR(0):      %+.3f μm\n", x_ir_0 * 1e6)
        @printf("  x_IR(end):    %+.3f μm\n", x_ir_end * 1e6)
        @printf("  x_IR range:   [%.1f, %.1f] μm\n", x_ir_min * 1e6, x_ir_max * 1e6)
    catch e
        @printf("  [Inner Race error: %s]\n", sprint(showerror, e))
    end

    try
        # ── Cage Speed ──
        println("\n--- Cage Speed (last 20% average) ---")
        i_ss = max(1, Int(round(0.8 * n_t)))
        cage_speeds = Float64[]
        for j in 1:Z
            θ_idx = ADORE.ball_pos_offset(j) + 2  # θ is 3rd component (0-indexed +2)
            Δθ = u_mat[end, θ_idx] - u_mat[i_ss, θ_idx]
            Δt = t_dim[end] - t_dim[i_ss]
            if Δt > 0
                push!(cage_speeds, Δθ / Δt * 30 / π)
            end
        end
        cage_theory = 0.5 * ω_ir * (1 - geom.d / geom.d_m * cos(geom.alpha_0)) * 30 / π
        @printf("  Cage speed:   %.1f ± %.1f RPM\n", mean(cage_speeds), std(cage_speeds))
        @printf("  Theory (OR):  %.1f RPM\n", cage_theory)
        @printf("  Ratio:        %.4f\n", mean(cage_speeds) / cage_theory)
    catch e
        @printf("  [Cage speed error: %s]\n", sprint(showerror, e))
    end

    try
        # ── Quaternion Integrity ──
        println("\n--- Quaternion Integrity ---")
        q_norms = Float64[]
        for j in 1:Z
            off = ADORE.ball_pos_offset(j) + 3
            qn = sqrt(sum(u_mat[end, off:off+3] .^ 2))
            push!(q_norms, qn)
        end
        @printf("  |q| range:    [%.10f, %.10f]\n", minimum(q_norms), maximum(q_norms))
        @printf("  max|q-1|:     %.2e\n", maximum(abs.(q_norms .- 1.0)))
    catch e
        @printf("  [Quaternion error: %s]\n", sprint(showerror, e))
    end

    try
        # ── Ball Spin Speed ──
        println("\n--- Ball 1 Spin (end) ---")
        bω_off = ADORE.ball_vel_offset(1, Z) + 3
        ω_b = [u_mat[end, bω_off+k] * scales.W for k in 0:2]
        @printf("  ω = [%.0f, %.0f, %.0f] rad/s\n", ω_b...)
        @printf("  |ω| = %.0f rad/s (%.0f RPM)\n", norm(ω_b), norm(ω_b) * 30 / π)
    catch e
        @printf("  [Ball spin error: %s]\n", sprint(showerror, e))
    end

    try
        # ── IR Velocity ──
        println("\n--- Inner Race Velocity (end) ---")
        ir_v_off = ADORE.ir_vel_offset(Z)
        v_ir = [u_mat[end, ir_v_off+k] * scales.V for k in 0:2]
        @printf("  v_IR = [%.3e, %.3e, %.3e] m/s\n", v_ir...)
        @printf("  |v_IR| = %.3e m/s\n", norm(v_ir))
    catch e
        @printf("  [IR velocity error: %s]\n", sprint(showerror, e))
    end

    println("\n--- Solution Statistics ---")
    @printf("  max|u★|:      %.3e\n", maximum(abs.(u_mat)))
    @printf("  State size:   %d DOFs\n", size(u_mat, 2))
end

println("\n" * "="^60)
println("  DONE at ", Dates.format(Dates.now(), "HH:MM:SS"))
println("="^60)
