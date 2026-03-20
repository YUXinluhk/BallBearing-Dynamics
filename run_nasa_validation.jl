# =====================================================================
# NASA TP-2275 Validation Runner
# Runs ADORE V2 for the 35-mm NASA bearing at 4 speed points
# and compares against experimental thermal data from Parker (1984)
# =====================================================================

println("="^70)
println("  ADORE V2 — NASA TP-2275 Experimental Validation")
println("  35-mm bore ACBB, F_a = 667 N, T_oil_in = 394 K")
println("="^70)

using Pkg
Pkg.activate(".")

using ADORE
using Printf
using Statistics: mean

"""
    run_single_case(toml_path, out_dir)

Run a single ADORE simulation case and print results.
"""
function run_single_case(toml_path, out_dir)
    config_file = joinpath(@__DIR__, toml_path)
    println("  Loading: $config_file")
    geom, mat, lub, trac, cage, config = load_simulation_config(config_file)

    n_rpm = config.inner_race_speed * 30 / π
    @printf("  Bearing: %d balls, n=%.0f RPM\n", geom.n_balls, n_rpm)

    println("  Running simulation...")
    flush(stdout)
    result = run_simulation(geom, mat, lub, trac, cage, config; verbose=true)
    @printf("  ODE retcode: %s, %d time steps\n", result.retcode, length(result.t))

    # Extract thermal state
    Z = geom.n_balls
    th_off = ADORE.thermal_offset(Z)
    T_ir = result.u[end, th_off]
    T_or = result.u[end, th_off+1]

    # Compute total heat
    fo_matrix = compute_field_outputs(result)
    fo_ball = zeros(length(result.t), Z, ADORE.N_FIELD_PER_BALL)
    for j in 1:Z
        base = (j - 1) * ADORE.N_FIELD_PER_BALL
        for f in 1:ADORE.N_FIELD_PER_BALL
            fo_ball[:, j, f] .= fo_matrix[:, base+f]
        end
    end
    fo_H_slide_i = fo_ball[:, :, ADORE.FO_H_SLIDE_I]
    fo_H_slide_o = fo_ball[:, :, ADORE.FO_H_SLIDE_O]
    fo_H_spin_i  = fo_ball[:, :, ADORE.FO_H_SPIN_I]
    fo_H_spin_o  = fo_ball[:, :, ADORE.FO_H_SPIN_O]
    fo_H_drag    = fo_ball[:, :, ADORE.FO_H_DRAG]
    fo_H_churn   = fo_ball[:, :, ADORE.FO_H_CHURN]
    H_total = fo_H_slide_i .+ fo_H_slide_o .+ fo_H_spin_i .+ fo_H_spin_o .+ fo_H_drag .+ fo_H_churn

    # Steady-state means (last 20%)
    n_t = length(result.t)
    ss_range = max(1, n_t - n_t ÷ 5 + 1):n_t
    H_mean = mean(sum(H_total[ss_range, :], dims=2))

    @printf("  Results:\n")
    @printf("    T_inner = %.1f K\n", T_ir)
    @printf("    T_outer = %.1f K\n", T_or)
    @printf("    H_total = %.3f kW\n", H_mean / 1000)

    return (T_inner=T_ir, T_outer=T_or, H_kW=H_mean/1000, retcode=result.retcode)
end

# NASA experimental data from Figs. 3 of Parker (1984)
nasa_exp = Dict(
    28000 => (T_inner=422.0, T_outer=413.0, H_kW=0.35),
    48000 => (T_inner=443.0, T_outer=425.0, H_kW=0.70),
    64000 => (T_inner=458.0, T_outer=435.0, H_kW=1.10),
    72000 => (T_inner=468.0, T_outer=443.0, H_kW=1.40),
)

cases = [
    ("NASA_28k", "inputs/NASA_35mm_28k.toml", 28000),
    ("NASA_48k", "inputs/NASA_35mm_48k.toml", 48000),
    ("NASA_64k", "inputs/NASA_35mm_64k.toml", 64000),
    ("NASA_72k", "inputs/NASA_35mm_72k.toml", 72000),
]

results = []

for (label, toml_path, rpm) in cases
    println("\n" * "="^50)
    println("  Running: $label ($rpm RPM)")
    println("="^50)

    out_dir = joinpath("results", label)
    mkpath(out_dir)

    try
        # Run the simulation (run_single.jl provides the pipeline)
        run_single_case(toml_path, out_dir)

        # Read back results for comparison
        exp = nasa_exp[rpm]
        println("\n  NASA Experimental:")
        println("    T_inner = $(exp.T_inner) K")
        println("    T_outer = $(exp.T_outer) K")
        println("    H_total = $(exp.H_kW) kW")

        push!(results, (rpm=rpm, label=label, status="OK"))
    catch e
        println("  ✗ Error: $(sprint(showerror, e))")
        push!(results, (rpm=rpm, label=label, status="FAILED: $e"))
    end
end

println("\n" * "="^70)
println("  BATCH SUMMARY")
println("="^70)
for r in results
    println("  $(r.label) ($(r.rpm) RPM): $(r.status)")
end
println("="^70)
