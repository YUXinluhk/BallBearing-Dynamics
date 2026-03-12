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

include("run_single.jl")  # reuse the existing single-case runner

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
