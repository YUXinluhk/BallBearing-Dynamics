# Diagnose error with full traceback
using ADORE
using Printf
using Statistics: mean

geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

config = SimulationConfig(
    t_end=0.005,     # short test: 5 ms
    dt_output=1e-3,
    inner_race_speed=10000.0 * pi / 30,
    F_axial=2000.0,
    t_ramp_end=5e-3,
    zeta=0.03,
    c_structural=20.0,
    mu_spin=0.06,
    integrator=IntegratorConfig(
        rtol=1e-5,
        atol=1e-7,
        h_max=1e-4,
        max_steps=1000,
    ),
)

println("Launching short simulation...")
try
    result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)
    println("SUCCESS! retcode = $(result.retcode)")
catch e
    println("\n=== ERROR ===")
    println(sprint(showerror, e, catch_backtrace()))
end
