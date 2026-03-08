#!/usr/bin/env julia
"""
Export Julia 30ms simulation results to CSV for comparison with Python.
"""

println("="^60)
println("  Julia — Export 30ms Results for Comparison")
println("="^60)
flush(stdout)

t_wall = time()

using ADORE
using Printf
using Statistics: mean, std
using LinearAlgebra: norm
using DelimitedFiles

geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

n_rpm = 10000.0
F_a = 2000.0
omega = n_rpm * π / 30

config = SimulationConfig(
    t_end=0.030,           # 30ms (same as Python)
    dt_output=50e-6,       # 50μs output (same as Python)
    inner_race_speed=omega,
    F_axial=F_a,
    t_ramp_end=0.005,      # 5ms cosine ramp (MATCHING Python)
    zeta=0.10,             # per-mode damping ratio (same as Python)
    c_structural=0.0,      # NOT USED — per-mode replaces this
    mu_spin=0.06,
    integrator=IntegratorConfig(
        rtol=1e-4,         # MATCHING Python
        atol=1e-7,         # same as Python
        h_max=1e-4,        # same as Python
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

# ── Extract time histories ──
# Inner race: x, y, z (positions, dimensional meters)
x_ir = u_mat[:, 1] .* scales.L
y_ir = u_mat[:, 2] .* scales.L
z_ir = u_mat[:, 3] .* scales.L

# Inner race velocity
ir_v_off = ADORE.ir_vel_offset(Z)
vx_ir = u_mat[:, ir_v_off] .* scales.V
vy_ir = u_mat[:, ir_v_off+1] .* scales.V
vz_ir = u_mat[:, ir_v_off+2] .* scales.V

# Cage speed (θ̇) → dimensional rad/s
cage_off = ADORE.cage_vel_offset(Z)
cage_omega = u_mat[:, cage_off+3] .* scales.W

# Ball orbital speeds — use velocity DOF θ̇★ directly (same method as cage ω)
ball_orbital = zeros(n_t, Z)
for j in 1:Z
    θ̇_idx = ADORE.ball_vel_offset(j, Z) + 2  # θ̇ velocity DOF
    ball_orbital[:, j] .= u_mat[:, θ̇_idx] .* scales.W  # nondim → rad/s
end
# Mean ball orbital speed
mean_ball_orbital = vec(mean(ball_orbital, dims=2))

# Ball 1 spin (|ω|)
ball1_spin = zeros(n_t)
for i in 1:n_t
    bω_off = ADORE.ball_vel_offset(1, Z) + 3
    ω_ball = [u_mat[i, bω_off+k] * scales.W for k in 0:2]
    ball1_spin[i] = norm(ω_ball)
end

# ── Export full state as numpy-compatible NPZ ──
# Write individual .npy files then Python reads them

function write_npy(path::String, arr::AbstractMatrix{Float64})
    open(path, "w") do io
        # numpy .npy format v1.0 header
        magic = UInt8[0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59]  # "\x93NUMPY"
        write(io, magic)
        write(io, UInt8(1))  # major version
        write(io, UInt8(0))  # minor version
        # Header dict (numpy uses Fortran order by default, we write C order)
        rows, cols = size(arr)
        header = "{'descr': '<f8', 'fortran_order': False, 'shape': ($rows, $cols), }        \n"
        hlen = UInt16(length(header))
        write(io, hlen)
        write(io, header)
        # Write data in C order (row-major)
        for i in 1:rows
            for j in 1:cols
                write(io, arr[i, j])
            end
        end
    end
end

function write_npy_vec(path::String, arr::AbstractVector{Float64})
    open(path, "w") do io
        magic = UInt8[0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59]
        write(io, magic)
        write(io, UInt8(1))
        write(io, UInt8(0))
        n = length(arr)
        header = "{'descr': '<f8', 'fortran_order': False, 'shape': ($n,), }                 \n"
        write(io, UInt16(length(header)))
        write(io, header)
        write(io, arr)
    end
end

function write_npy_scalar(path::String, val::Float64)
    write_npy_vec(path, [val])
end

# Save to npy directory for Python postprocessing
npy_dir = "julia_npy"
mkpath(npy_dir)

write_npy_vec(joinpath(npy_dir, "t.npy"), t_dim)
write_npy(joinpath(npy_dir, "u_mat.npy"), u_mat)
write_npy_scalar(joinpath(npy_dir, "L.npy"), scales.L)
write_npy_scalar(joinpath(npy_dir, "V.npy"), scales.V)
write_npy_scalar(joinpath(npy_dir, "W.npy"), scales.W)
write_npy_scalar(joinpath(npy_dir, "Q.npy"), scales.Q)

@printf("\n  Exported %d rows to %s/ (npy format)\n", n_t, npy_dir)

# Also keep CSV for quick inspection
outfile = "julia_30ms_results.csv"
header = "t_s,x_ir_m,y_ir_m,z_ir_m,vx_ir_ms,cage_omega_rads,mean_ball_orbital_rads,ball1_spin_rads"
out_csv = hcat(
    t_dim,
    u_mat[:, 1] .* scales.L,
    u_mat[:, 2] .* scales.L,
    u_mat[:, 3] .* scales.L,
    u_mat[:, ir_v_off] .* scales.V,
    cage_omega,
    mean_ball_orbital,
    ball1_spin,
)
open(outfile, "w") do io
    println(io, header)
    writedlm(io, out_csv, ',')
end
@printf("  Exported %d rows to %s (CSV)\n", n_t, outfile)

# Theoretical cage speed
cage_theory_rads = 0.5 * omega * (1 - geom.d / geom.d_m * cos(geom.alpha_0))
@printf("  Theory cage speed: %.1f rad/s = %.1f RPM\n", cage_theory_rads, cage_theory_rads * 30 / π)
@printf("  Final cage speed:  %.1f rad/s = %.1f RPM (%.1f%%)\n",
    cage_omega[end], cage_omega[end] * 30 / π, cage_omega[end] / cage_theory_rads * 100)

println("\n" * "="^60)
println("  EXPORT COMPLETE")
println("="^60)

