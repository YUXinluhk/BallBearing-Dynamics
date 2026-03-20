# =====================================================================
# Dynamics/driver.jl — Simulation orchestrator
#
# Port of dynamics.py::BearingDynamics.run.
# Wires together: config → Hertz precompute → scales → params → IC →
# ODE problem → solve → result.
# =====================================================================

using OrdinaryDiffEq: ODEProblem, ODEFunction, solve, Tsit5, FBDF
using ADTypes: AutoForwardDiff  # 极致性能：全部奇异点被消灭，重新启用精确解析 Jacobian !
using DiffEqCallbacks
using Statistics: mean, std
using Printf

"""
    SimResult — Output of transient dynamic simulation.
"""
struct SimResult
    t::Vector{Float64}          # time samples [s] (dimensional)
    u::Matrix{Float64}          # state matrix (n_samples × n_state)
    scales::Scales
    params::ODEParams
    ca_axes::Any                # ComponentArray axes for state vector reconstruction
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
            bp = u.ball[j].pos
            s = bp.q0
            v1 = bp.q1
            v2 = bp.q2
            v3 = bp.q3
            n = sqrt(s^2 + v1^2 + v2^2 + v3^2)
            if n > 1e-15
                bp.q0 = s / n
                bp.q1 = v1 / n
                bp.q2 = v2 / n
                bp.q3 = v3 / n
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

    params_tuple = (params,)  # wrap ODEParams for kernel

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

    sol = solve(prob, FBDF(; autodiff=AutoForwardDiff());
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

    return SimResult(t_dim, u_mat, scales, params, getaxes(u0_star), sol.retcode)
end

"""
    run_convergence_scan(geom, mat, lub, trac, cage, config;
        h_max_values=[2e-6, 1e-6, 5e-7, 2e-7, 1e-7],
        symmetry_threshold=1e-3) → DataFrame-like output

Run the same simulation with different `h_max` values and report symmetry metrics.
Returns a vector of NamedTuples with convergence diagnostics.

Used for TC-PHY-006 numerical robustness validation.
"""
function run_convergence_scan(geom::BearingGeometry, mat::MaterialParams,
    lub::LubricantParams, trac::TractionParams,
    cage::CageGeometry, config::SimulationConfig;
    h_max_values::Vector{Float64}=[2e-6, 1e-6, 5e-7, 2e-7, 1e-7],
    symmetry_threshold::Float64=1e-3)

    Z = geom.n_balls
    results = []

    for hm in h_max_values
        ic_mod = IntegratorConfig(
            rtol=config.integrator.rtol,
            atol=config.integrator.atol,
            h_max=hm,
            max_steps=config.integrator.max_steps,
            eps_contact=config.integrator.eps_contact
        )
        config_mod = SimulationConfig(
            t_end=config.t_end,
            dt_output=config.dt_output,
            inner_race_speed=config.inner_race_speed,
            outer_race_speed=config.outer_race_speed,
            F_axial=config.F_axial,
            F_radial=config.F_radial,
            t_ramp_end=config.t_ramp_end,
            mu_spin=config.mu_spin,
            c_structural=config.c_structural,
            zeta=config.zeta,
            alpha2_relax_time=config.alpha2_relax_time,
            delta_r_thermal=config.delta_r_thermal,
            integrator=ic_mod,
            churning=config.churning,
            thermal=config.thermal
        )

        @printf("  h_max = %.1e ... ", hm)
        flush(stdout)

        sim = run_simulation(geom, mat, lub, trac, cage, config_mod)

        # Compute per-ball Q_i symmetry at each saved time
        n_t = length(sim.t)
        first_qi_gt = NaN
        first_qo_gt = NaN
        peak_qi_std = 0.0
        peak_qo_std = 0.0
        peak_qi_time = 0.0
        peak_qo_time = 0.0

        for k in 1:n_t
            u_k = sim.u[k, :]
            fo = compute_field_outputs(u_k, sim.params, Z)
            qi_vals = [fo[j].Q_i_dim for j in 1:Z]
            qo_vals = [fo[j].Q_o_dim for j in 1:Z]
            qi_std = std(qi_vals)
            qo_std = std(qo_vals)

            if qi_std > symmetry_threshold && isnan(first_qi_gt)
                first_qi_gt = sim.t[k]
            end
            if qo_std > symmetry_threshold && isnan(first_qo_gt)
                first_qo_gt = sim.t[k]
            end
            if qi_std > peak_qi_std
                peak_qi_std = qi_std
                peak_qi_time = sim.t[k]
            end
            if qo_std > peak_qo_std
                peak_qo_std = qo_std
                peak_qo_time = sim.t[k]
            end
        end

        @printf("PEAK_QI_STD=%.4e  FIRST_GT=%.4e s\n", peak_qi_std,
            isnan(first_qi_gt) ? NaN : first_qi_gt)

        push!(results, (
            h_max=hm,
            retcode=sim.retcode,
            n_steps=n_t,
            first_qi_gt=first_qi_gt,
            first_qo_gt=first_qo_gt,
            peak_qi_std=peak_qi_std,
            peak_qo_std=peak_qo_std,
            peak_qi_time=peak_qi_time,
            peak_qo_time=peak_qo_time,
        ))
    end

    # Print summary table
    println("\n  ── Convergence Scan Summary ──")
    @printf("  %-10s  %-14s  %-14s  %-12s  %-12s\n",
        "h_max", "FIRST_QI_GT", "FIRST_QO_GT", "PEAK_QI_STD", "PEAK_QO_STD")
    for r in results
        @printf("  %-10.1e  %-14s  %-14s  %-12.4e  %-12.4e\n",
            r.h_max,
            isnan(r.first_qi_gt) ? "NaN" : @sprintf("%.4e s", r.first_qi_gt),
            isnan(r.first_qo_gt) ? "NaN" : @sprintf("%.4e s", r.first_qo_gt),
            r.peak_qi_std, r.peak_qo_std)
    end

    return results
end
