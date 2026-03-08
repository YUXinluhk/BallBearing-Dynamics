using ADORE
using Printf

geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

n_rpm = 10000.0
F_a = 2000.0
omega = n_rpm * pi / 30

config = SimulationConfig(
    t_end=0.005,
    dt_output=100e-6,
    inner_race_speed=omega,
    F_axial=F_a,
    t_ramp_end=0.005,
    zeta=0.10,
    c_structural=0.0,
    mu_spin=0.06,
    integrator=IntegratorConfig(
        rtol=1e-4, atol=1e-7, h_max=1e-4, max_steps=10_000_000
    ),
)

result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)

Z = geom.n_balls
scales = result.scales
t_dim = result.t
u_mat = result.u

cage_off = ADORE.cage_vel_offset(Z)
cage_theory = 0.5 * omega * (1 - geom.d * cos(geom.alpha_0) / geom.d_m)

for i in [1, div(length(t_dim), 4), div(length(t_dim), 2), length(t_dim)]
    cage_w = u_mat[i, cage_off+3] * scales.W
    ball_vx = u_mat[i, ADORE.ball_vel_offset(1, Z)] * scales.V
    x_ir = u_mat[i, 1] * scales.L
    @printf("t=%.4f ms  cage=%.1f RPM (%.1f%%)  ball_vx=%.2e m/s  x_ir=%.2e m\n",
        t_dim[i] * 1e3, cage_w * 30 / pi, cage_w / cage_theory * 100, ball_vx, x_ir)
end
