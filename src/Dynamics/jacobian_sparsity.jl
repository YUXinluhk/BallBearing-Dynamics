# =====================================================================
# Dynamics/jacobian_sparsity.jl — Analytical Jacobian sparsity pattern
#
# Derived from ComponentArray layout (not flat offset calculators).
# Ball j only couples to its own states + inner race + cage + thermal.
# =====================================================================

"""
    build_jacobian_sparsity(Z::Int) → SparseMatrixCSC{Float64}

Build the analytical Jacobian sparsity pattern for the bearing ODE.

Uses a dummy ComponentArray to determine actual flat indices, ensuring
the pattern matches the ComponentArray memory layout used by the kernel.
"""
function build_jacobian_sparsity(Z::Int)
    # Build a dummy CA with same structure as init_state creates
    ball_states = [(pos=(x=0.0, r=0.0, θ=0.0, q0=0.0, q1=0.0, q2=0.0, q3=0.0),
                    vel=(x=0.0, r=0.0, θ=0.0, ωx=0.0, ωy=0.0, ωz=0.0)) for _ in 1:Z]
    dummy_nt = (
        ir = (pos=(x=0.0, y=0.0, z=0.0, γy=0.0, γz=0.0), vel=(x=0.0, y=0.0, z=0.0, γy=0.0, γz=0.0)),
        cage = (pos=(x=0.0, y=0.0, z=0.0, θ=0.0), vel=(x=0.0, y=0.0, z=0.0, θ=0.0)),
        ball = ball_states,
        thermal = (T_i=0.0, T_o=0.0, T_b=0.0, T_oil=0.0),
        heat = (ir=0.0, or_=0.0, b=0.0, oil=0.0, Q_amb=0.0, Cp=0.0)
    )
    ca = ComponentArray(dummy_nt)
    N = length(ca)

    # Helper: get flat index range of a ComponentArray view
    function ca_indices(view_expr)
        # Mark the view with sequential indices, read them back
        test = ComponentArray(collect(1.0:N), getaxes(ca)...)
        v = view_expr(test)
        return Int.(collect(v))
    end

    rows = Int[]
    cols = Int[]

    function mark_block!(r_range, c_range)
        for r in r_range, c in c_range
            push!(rows, r)
            push!(cols, c)
        end
    end

    # ── Index ranges from ComponentArray ──
    ir_p  = ca_indices(u -> u.ir.pos)
    ir_v  = ca_indices(u -> u.ir.vel)
    cage_p = ca_indices(u -> u.cage.pos)
    cage_v = ca_indices(u -> u.cage.vel)
    ball_p = [ca_indices(u -> u.ball[j].pos) for j in 1:Z]
    ball_v = [ca_indices(u -> u.ball[j].vel) for j in 1:Z]
    th    = ca_indices(u -> u.thermal)
    ha    = ca_indices(u -> u.heat)

    # ── 1. Kinematic: position derivatives = velocities ──
    mark_block!(ir_p, ir_v)
    mark_block!(cage_p, cage_v)
    for j in 1:Z
        mark_block!(ball_p[j], ball_v[j])
        mark_block!(ball_p[j], ball_p[j])   # quaternion self-coupling
    end

    # ── 2. Force coupling: velocity derivatives ──
    for j in 1:Z
        mark_block!(ball_v[j], ball_p[j])    # ball pos → ball acc
        mark_block!(ball_v[j], ball_v[j])    # ball vel → ball acc (damping)
        mark_block!(ball_v[j], ir_p)         # IR pos → ball acc
        mark_block!(ball_v[j], ir_v)         # IR vel → ball acc (surface vel)
        mark_block!(ball_v[j], cage_p)       # cage pos → ball acc (pocket)
        mark_block!(ball_v[j], cage_v)       # cage vel → ball acc (pocket)

        # Newton's 3rd law
        mark_block!(ir_v, ball_p[j])
        mark_block!(ir_v, ball_v[j])
        mark_block!(cage_v, ball_p[j])
        mark_block!(cage_v, ball_v[j])
    end

    # Self-coupling
    mark_block!(ir_v, ir_p)
    mark_block!(ir_v, ir_v)
    mark_block!(cage_v, cage_p)
    mark_block!(cage_v, cage_v)

    # ── 3. Thermal coupling ──
    mark_block!(th, th)
    mark_block!(ha, ha)

    # Temperature → friction → velocity
    mark_block!(ir_v, th)
    mark_block!(cage_v, th)
    for j in 1:Z
        mark_block!(ball_v[j], th)
    end

    # Cage pilot ↔ IR
    mark_block!(ir_v, cage_p)
    mark_block!(ir_v, cage_v)
    mark_block!(cage_v, ir_p)
    mark_block!(cage_v, ir_v)

    # Build sparse matrix
    vals = ones(Float64, length(rows))
    S = sparse(rows, cols, vals, N, N)

    return S
end
