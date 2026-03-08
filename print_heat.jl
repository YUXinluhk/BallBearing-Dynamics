using ADORE
using Printf
using Statistics

geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage = cage_from_bearing(geom)

config = SimulationConfig(
    t_end=0.005, dt_output=50e-6,
    inner_race_speed=10000.0 * pi / 30, F_axial=2000.0,
    t_ramp_end=0.005, zeta=0.10, c_structural=0.0, mu_spin=0.06
)

res = run_simulation(geom, mat, lub, trac, cage, config, verbose=false)
fo = compute_field_outputs(res)

Z = geom.n_balls
n_t = length(res.t)
n_steady = max(1, n_t ÷ 5)
ss_range = max(1, n_t - n_steady + 1):n_t

fo_ball = zeros(n_t, Z, ADORE.N_FIELD_PER_BALL)
for j in 1:Z
    base = (j - 1) * ADORE.N_FIELD_PER_BALL
    for f in 1:ADORE.N_FIELD_PER_BALL
        fo_ball[:, j, f] .= fo[:, base+f]
    end
end

labels = ["Inner Sliding", "Outer Sliding", "Inner Spin", "Outer Spin", "Churning Drag", "Churning Rot."]
idxs = [ADORE.FO_H_SLIDE_I, ADORE.FO_H_SLIDE_O, ADORE.FO_H_SPIN_I, ADORE.FO_H_SPIN_O, ADORE.FO_H_DRAG, ADORE.FO_H_CHURN]

println("\n--- Heat Breakdown (W) ---")
total_heat = 0.0
for (lbl, idx) in zip(labels, idxs)
    global total_heat
    val = mean(sum(fo_ball[ss_range, :, idx], dims=2))
    @printf("%-20s: %8.3f\n", lbl, val)
    total_heat += val
end
@printf("%-20s: %8.3f\n", "TOTAL", total_heat)
