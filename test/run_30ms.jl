#!/usr/bin/env julia
"""
Julia — 30ms Dynamic Simulation
"""

println("="^60)
println("  Julia — 30ms Dynamic Simulation")
println("="^60)
flush(stdout)

t_wall = time()

using ADORE
using Printf
using Statistics: mean, std
using LinearAlgebra: norm

geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

n_rpm = 10000.0
F_a = 2000.0
omega = n_rpm * π / 30

config = SimulationConfig(
    t_end=0.030,           # 30ms
    dt_output=50e-6,       # 50μs output → 600 points
    inner_race_speed=omega,
    F_axial=F_a,
    t_ramp_end=0.0,        # NO ramp — ICs from QS at full speed
    zeta=0.10,             # per-mode damping ratio (matching Python)
    c_structural=0.0,      # NOT USED — per-mode replaces this
    mu_spin=0.06,
    integrator=IntegratorConfig(
        rtol=1e-5,
        atol=1e-7,
        h_max=1e-4,
        max_steps=10_000_000,
    ),
)

@printf("  Bearing: 7210B (%d balls)\n", geom.n_balls)
@printf("  Load: Fa=%.0f N, n=%.0f RPM\n", F_a, n_rpm)
@printf("  Duration: %.0f ms\n", config.t_end * 1000)
flush(stdout)

println("\nRunning...")
flush(stdout)

result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)

elapsed = time() - t_wall
@printf("\n  Wall time: %.1f s\n", elapsed)
@printf("  ODE retcode: %s\n", result.retcode)
@printf("  Time points: %d\n", length(result.t))

t_dim = result.t
u_mat = result.u
Z = geom.n_balls
scales = result.scales
n_t = length(t_dim)

has_nan = any(isnan.(u_mat))
@printf("  NaN: %s\n", has_nan)

if !has_nan && n_t > 5
    # Inner race axial position
    x_ir = u_mat[:, 1] .* scales.L
    @printf("\n--- Inner Race Axial Position ---\n")
    @printf("  x_IR(0)   = %+.3f μm\n", x_ir[1] * 1e6)
    @printf("  x_IR(end) = %+.3f μm\n", x_ir[end] * 1e6)
    @printf("  x_IR range: [%.1f, %.1f] μm\n", minimum(x_ir) * 1e6, maximum(x_ir) * 1e6)

    # Last 10ms statistics (settled region)
    last_10ms_idx = findall(t_dim .>= 0.020)
    if length(last_10ms_idx) > 2
        x_ir_last = x_ir[last_10ms_idx]
        @printf("  x_IR(20-30ms) mean: %.2f μm, std: %.3f μm\n",
            mean(x_ir_last) * 1e6, std(x_ir_last) * 1e6)
    end

    # Ball orbital speeds (cage-like speed)
    @printf("\n--- Ball Orbital Speeds (t_end) ---\n")
    ball_rpm = Float64[]
    for j in 1:Z
        θ_idx = ADORE.ball_pos_offset(j) + 2  # θ position
        if n_t > 2
            Δθ = u_mat[end, θ_idx] - u_mat[end-1, θ_idx]
            Δt = t_dim[end] - t_dim[end-1]
            if Δt > 0
                ωj = Δθ / Δt * 30 / π
                push!(ball_rpm, ωj)
            end
        end
    end
    cage_theory = 0.5 * omega * (1 - geom.d / geom.d_m * cos(geom.alpha_0)) * 30 / π
    if !isempty(ball_rpm)
        @printf("  Mean ball θ̇: %.1f RPM\n", mean(ball_rpm))
        @printf("  Theory cage: %.1f RPM\n", cage_theory)
        @printf("  Ratio: %.4f\n", mean(ball_rpm) / cage_theory)
    end

    # Ball 1 spin
    @printf("\n--- Ball 1 Spin (t_end) ---\n")
    bω_off = ADORE.ball_vel_offset(1, Z) + 3
    ω_ball = [u_mat[end, bω_off+k] * scales.W for k in 0:2]
    @printf("  ω = [%.0f, %.0f, %.0f] rad/s\n", ω_ball...)
    @printf("  |ω| = %.0f rad/s (%.0f RPM)\n", norm(ω_ball), norm(ω_ball) * 30 / π)

    # Inner race velocity
    @printf("\n--- Inner Race Velocity (t_end) ---\n")
    ir_v_off = ADORE.ir_vel_offset(Z)
    v_ir = [u_mat[end, ir_v_off+k] * scales.V for k in 0:2]
    @printf("  v_IR = [%.3e, %.3e, %.3e] m/s\n", v_ir...)
    @printf("  |v_IR| = %.3e m/s\n", norm(v_ir))

    # Cage speed
    @printf("\n--- Cage Speed (t_end) ---\n")
    cage_off = ADORE.cage_vel_offset(Z)
    cage_ωθ = u_mat[end, cage_off+3] * scales.W  # cage θ̇
    cage_rpm = cage_ωθ * 30 / π
    @printf("  ω_cage = %.1f rad/s = %.1f RPM\n", cage_ωθ, cage_rpm)
    @printf("  Theory: %.1f RPM\n", cage_theory)
    @printf("  Ratio: %.4f\n", cage_rpm / cage_theory)

    # Quaternion check
    @printf("\n--- Quaternion Integrity ---\n")
    q_norms = Float64[]
    for j in 1:Z
        off = ADORE.ball_pos_offset(j) + 3
        qn = sqrt(sum(u_mat[end, off:off+3] .^ 2))
        push!(q_norms, qn)
    end
    @printf("  |q| range: [%.10f, %.10f]\n", minimum(q_norms), maximum(q_norms))
    @printf("  max|q-1|: %.2e\n", maximum(abs.(q_norms .- 1.0)))
end

println("\n" * "="^60)
println("  JULIA 30ms TEST COMPLETE")
println("="^60)
