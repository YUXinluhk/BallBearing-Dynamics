# Pinpoint ForwardDiff typeassert error
using ADORE
using ForwardDiff: Dual, Tag

# Build same setup as ODE solver
geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

config = SimulationConfig(
    t_end=0.005,
    dt_output=1e-3,
    inner_race_speed=10000.0 * pi / 30,
    F_axial=2000.0,
    t_ramp_end=5e-3,
    zeta=0.03,
    c_structural=20.0,
    mu_spin=0.06,
    integrator=IntegratorConfig(rtol=1e-5, atol=1e-7, h_max=1e-4, max_steps=1000),
)

h_inner, h_outer = ADORE.create_bearing_hertz(geom, mat)
scales = ADORE.Scales(geom, h_inner)
qs = ADORE.solve_quasi_static(geom, mat; F_a=config.F_axial, F_rz=config.F_radial,
    n_rpm=config.inner_race_speed * 30 / pi, verbose=false)
params = ADORE.build_params(geom, mat, lub, trac, cage, config, scales, h_inner, h_outer; qs=qs)
u0 = ADORE.init_state(geom, qs, config)
u0_star = ADORE.nondim_state(u0, scales, geom.n_balls)

# Try calling ode_rhs! with Dual numbers
N = length(u0_star)
println("State vector length: $N")

# Create Dual type input
tag = Tag{Nothing,Float64}
u_dual = [Dual{tag}(u0_star[i], ntuple(j -> j == i ? 1.0 : 0.0, 1)...) for i in 1:N]
du_dual = similar(u_dual)

params_tuple = (params,)
t_star = 0.0

println("Calling ode_rhs! with Float64...")
du_f64 = similar(u0_star)
ADORE.ode_rhs!(du_f64, u0_star, params_tuple, t_star)
println("Float64 call succeeded!")

println("Calling ode_rhs! with Dual numbers (1 partial)...")
try
    ADORE.ode_rhs!(du_dual, u_dual, params_tuple, t_star)
    println("Dual call succeeded!")
catch e
    println("\n=== DUAL CALL FAILED ===")
    println(sprint(showerror, e, catch_backtrace()))
end
