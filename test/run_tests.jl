using ADORE, Test, Printf, Statistics

println("=== STAGE 1: QS ===")
geom, mat = bearing_7210B()
h_inner, h_outer = create_bearing_hertz(geom, mat)

qs = solve_quasi_static(geom, mat; F_a=2000.0, n_rpm=10000.0, verbose=true)
@printf("  Converged: %s, α_i=%.2f°, Q_i=%.1f N\n", qs.converged, mean(qs.alpha_inner), mean(qs.Q_inner))

println("\n=== STAGE 2: Scales & Params ===")
scales = Scales(geom, h_inner)
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)
config = SimulationConfig(t_end=0.001, dt_output=1e-5,
    inner_race_speed=10000.0 * π / 30, F_axial=2000.0)
params = build_params(geom, mat, lub, trac, cage, config, scales, h_inner, h_outer)
println("  Params OK, len=", length(params))

println("\n=== STAGE 3: State ===")
Z = geom.n_balls
try
    u0 = init_state(geom, qs, config)
    println("  State OK, len=", length(u0))
catch e
    println("  init_state FAILED: ", sprint(showerror, e))
end

println("\n=== STAGE 4: ODE RHS ===")
try
    u0 = init_state(geom, qs, config)
    u0_star = nondim_state(u0, scales, Z)
    du = similar(u0_star)
    ode_rhs!(du, u0_star, (params,), 0.0)
    @printf("  max|du| = %.3e, NaN=%s\n", maximum(abs.(du)), any(isnan.(du)))
catch e
    println("  ODE RHS FAILED: ", sprint(showerror, e))
end

println("\n=== DONE ===")
