# =====================================================================
# Dynamics/driver.jl — Simulation orchestrator
#
# Port of dynamics.py::BearingDynamics.run.
# Wires together: config → Hertz precompute → scales → params → IC →
# ODE problem → solve → result.
# =====================================================================

using OrdinaryDiffEq: ODEProblem, ODEFunction, solve, Tsit5, FBDF
using ADTypes: AutoFiniteDiff  # ForwardDiff blocked by SciMLBase.UJacobianWrapper AbstractFloat assertion; kernel is AD-ready
using DiffEqCallbacks
using Statistics: mean
using Printf

"""
    SimResult — Output of transient dynamic simulation.
"""
struct SimResult
    t::Vector{Float64}          # time samples [s] (dimensional)
    u::Matrix{Float64}          # state matrix (n_state × n_samples)
    scales::Scales
    params::Vector{Float64}
    retcode::Any
end

"""
    quaternion_renormalize_callback(Z)

Create a `DiscreteCallback` that renormalizes all ball quaternions every step.
"""
function quaternion_renormalize_callback(Z::Int)
    function affect!(integrator)
        u = integrator.u
        @inbounds for j in 1:Z
            off = ball_pos_offset(j) + 3
            s = u[off]
            v1 = u[off+1]
            v2 = u[off+2]
            v3 = u[off+3]
            n = sqrt(s^2 + v1^2 + v2^2 + v3^2)
            if n > 1e-15
                u[off] = s / n
                u[off+1] = v1 / n
                u[off+2] = v2 / n
                u[off+3] = v3 / n
            end
        end
    end
    DiscreteCallback((u, t, integrator) -> true, affect!)
end

"""
    run_simulation(geom, mat, lub, trac, cage, config; verbose=false) → SimResult

Run full 6-DOF transient bearing dynamics simulation.

Steps:
1. Hertz precomputation
2. Nondimensionalization scales
3. Quasi-static initial conditions
4. Build ODE problem
5. Integrate with OrdinaryDiffEq.jl
6. Return dimensional results
"""
function run_simulation(geom::BearingGeometry, mat::MaterialParams,
    lub::LubricantParams, trac::TractionParams,
    cage::CageGeometry, config::SimulationConfig;
    verbose::Bool=false)

    Z = geom.n_balls

    # 1. Hertz precompute
    h_inner, h_outer = create_bearing_hertz(geom, mat)
    verbose && println("Hertz inner Υ = $(h_inner.Upsilon) N/m^1.5")
    verbose && println("Hertz outer Υ = $(h_outer.Upsilon) N/m^1.5")

    # 2. Scales
    scales = Scales(geom, h_inner)
    verbose && println("Scales: Q=$(scales.Q) N, L=$(scales.L) m, T=$(scales.T) s")

    # 3. Quasi-static IC
    verbose && println("Solving quasi-static equilibrium...")
    qs = solve_quasi_static(geom, mat;
        F_a=config.F_axial,
        F_rz=config.F_radial,
        n_rpm=config.inner_race_speed * 30 / π,
        verbose=verbose)
    verbose && println("  QS converged: $(qs.converged)")
    verbose && println("  QS α_inner mean: $(mean(qs.alpha_inner))°")

    # 4. Build parameters
    params = build_params(geom, mat, lub, trac, cage, config, scales, h_inner, h_outer; qs=qs)

    # 5. Initial state
    u0 = init_state(geom, qs, config, lub)
    u0_star = nondim_state(u0, scales, Z)

    # 6. ODE Problem
    t_end_star = nondim_time(scales, config.t_end)
    tspan = (0.0, t_end_star)

    # Output times
    dt_star = nondim_time(scales, config.dt_output)
    saveat = dt_star

    params_tuple = (params,)  # wrap for kernel

    # Build sparse Jacobian prototype for efficient finite differences
    # jac_prototype goes on ODEFunction, not the algorithm constructor
    jac_sparsity = build_jacobian_sparsity(Z)
    N_dof = n_state(Z)
    verbose && @printf("  Jacobian: %d DOFs, %d/%d nonzeros (%.1f%% sparse)\n",
        N_dof, nnz(jac_sparsity), N_dof^2, 100.0 * (1 - nnz(jac_sparsity) / N_dof^2))

    f = ODEFunction(ode_rhs!; jac_prototype=float.(jac_sparsity))
    prob = ODEProblem(f, u0_star, tspan, params_tuple)

    # Callbacks: progress only (Baumgarte replaces quat callback)
    # Progress logger callback
    t_last_print = Ref(time())
    function print_progress(integrator)
        if time() - t_last_print[] > 2.0
            t_dim = integrator.t * scales.T
            perc = t_dim / config.t_end * 100
            @printf("  [ODE] Integrator at t = %7.5e s  (%.2f %%)\n", t_dim, perc)
            flush(stdout)
            t_last_print[] = time()
        end
    end
    # Throttle progress to every 1e-4 s (nondim) instead of every step
    dt_prog = nondim_time(scales, 1e-4)
    cb_prog = PeriodicCallback(print_progress, dt_prog)

    # 7. Integrate
    ic = config.integrator
    verbose && println("Integrating ($(config.t_end*1e3) ms)...")

    sol = solve(prob, FBDF(; autodiff=AutoFiniteDiff());
        reltol=ic.rtol, abstol=ic.atol,
        dtmax=nondim_time(scales, ic.h_max),
        maxiters=ic.max_steps,
        saveat=saveat,
        callback=cb_prog,
        progress=verbose)

    verbose && println("  ODE retcode: $(sol.retcode)")
    verbose && println("  $(length(sol.t)) time points saved")

    # 8. Package results
    t_dim = dim_time.(Ref(scales), sol.t)
    u_mat = hcat(sol.u...)'  # (n_time × n_state)

    return SimResult(t_dim, u_mat, scales, params, sol.retcode)
end
