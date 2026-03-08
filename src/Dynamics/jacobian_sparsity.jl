# =====================================================================
# Dynamics/jacobian_sparsity.jl — Analytical Jacobian sparsity pattern
#
# The bearing ODE has block-sparse structure: ball j only couples to
# its own states + inner race + cage (NOT to other balls).
# Declaring this pattern lets FiniteDiff color the Jacobian with
# ~40 groups instead of 226 columns → 5-10× fewer RHS evaluations.
# =====================================================================
#
# Requires: SparseArrays (imported at module level)

"""
    build_jacobian_sparsity(Z::Int) → SparseMatrixCSC{Float64}

Build the analytical Jacobian sparsity pattern for the bearing ODE.

The state vector layout is:
  Position: IR(5) + Ball(7)×Z + Cage(4)
  Velocity: IR(5) + Ball(6)×Z + Cage(4)
  Total = 18 + 13Z

Coupling structure:
  - Inner race ↔ all balls (via contact forces)
  - Cage ↔ all balls (via pocket forces)
  - Ball j ↔ only itself + IR + Cage (NOT other balls)
  - Position derivatives depend on velocities (kinematic)
  - Velocity derivatives depend on positions (forces) + velocities (damping)
"""
function build_jacobian_sparsity(Z::Int)
    N = n_state(Z)
    N_p = n_pos_dofs(Z)

    # Collect (row, col) pairs for non-zero entries
    rows = Int[]
    cols = Int[]

    function mark_block!(r_range, c_range)
        for r in r_range, c in c_range
            push!(rows, r)
            push!(cols, c)
        end
    end

    # ── Index ranges ──
    # Position ranges
    ir_p = ir_pos_offset():ir_pos_offset()+N_IR_POS-1
    cage_p = cage_pos_offset(Z):cage_pos_offset(Z)+N_CAGE_POS-1
    ball_p = [ball_pos_offset(j):ball_pos_offset(j)+N_BALL_POS-1 for j in 1:Z]

    # Velocity ranges
    ir_v = ir_vel_offset(Z):ir_vel_offset(Z)+N_IR_VEL-1
    cage_v = cage_vel_offset(Z):cage_vel_offset(Z)+N_CAGE_VEL-1
    ball_v = [ball_vel_offset(j, Z):ball_vel_offset(j, Z)+N_BALL_VEL-1 for j in 1:Z]

    # ── 1. Kinematic: position derivatives = velocities ──
    # d(IR_pos)/dt depends on IR_vel
    mark_block!(ir_p, ir_v)
    # d(cage_pos)/dt depends on cage_vel
    mark_block!(cage_p, cage_v)
    # d(ball_pos)/dt depends on ball_vel (+ quaternion depends on ω)
    for j in 1:Z
        mark_block!(ball_p[j], ball_v[j])
        # Quaternion also depends on own quaternion (q̇ = ½ q ⊗ ω)
        mark_block!(ball_p[j], ball_p[j])
    end

    # ── 2. Force coupling: velocity derivatives ──
    for j in 1:Z
        # Ball acceleration depends on:
        #   - own position (contact geometry, traction, drag)
        #   - own velocity (damping, Coriolis)
        #   - IR position (groove center → contact geometry)
        #   - IR velocity (surface velocity → slide speed)
        #   - cage position (pocket interaction)
        #   - cage velocity (pocket interaction)
        mark_block!(ball_v[j], ball_p[j])    # ball pos → ball acc
        mark_block!(ball_v[j], ball_v[j])    # ball vel → ball acc (damping)
        mark_block!(ball_v[j], ir_p)         # IR pos → ball acc
        mark_block!(ball_v[j], ir_v)         # IR vel → ball acc (surface vel)
        mark_block!(ball_v[j], cage_p)       # cage pos → ball acc (pocket)
        mark_block!(ball_v[j], cage_v)       # cage vel → ball acc (pocket)

        # IR acceleration depends on ball j (Newton's 3rd law):
        mark_block!(ir_v, ball_p[j])         # ball pos → IR acc
        mark_block!(ir_v, ball_v[j])         # ball vel → IR acc

        # Cage acceleration depends on ball j (pocket force):
        mark_block!(cage_v, ball_p[j])       # ball pos → cage acc
        mark_block!(cage_v, ball_v[j])       # ball vel → cage acc
    end

    # IR acceleration depends on own state (damping)
    mark_block!(ir_v, ir_p)
    mark_block!(ir_v, ir_v)

    # Cage acceleration depends on own state (pilot, damping)
    mark_block!(cage_v, cage_p)
    mark_block!(cage_v, cage_v)

    # Build sparse matrix
    vals = ones(Float64, length(rows))
    S = sparse(rows, cols, vals, N, N)

    return S
end
