# Quick test to find the right imports
println("Testing imports...")

try
    using OrdinaryDiffEq: FBDF
    println("FBDF available: ", FBDF)
catch e
    println("FBDF not available: ", e)
end

try
    using ADTypes: AutoFiniteDiff
    println("ADTypes.AutoFiniteDiff: ", AutoFiniteDiff)
catch e
    println("ADTypes.AutoFiniteDiff not available: ", e)
end

try
    using OrdinaryDiffEq: QNDF
    println("QNDF available: ", QNDF)
catch e
    println("QNDF not available: ", e)
end

# List what's available
using OrdinaryDiffEq
println("\nHas FBDF: ", isdefined(OrdinaryDiffEq, :FBDF))
println("Has QNDF: ", isdefined(OrdinaryDiffEq, :QNDF))
println("Has KenCarp4: ", isdefined(OrdinaryDiffEq, :KenCarp4))
println("Has ImplicitEuler: ", isdefined(OrdinaryDiffEq, :ImplicitEuler))
