using Test
using ADORE

@testset "ADORE.jl" begin
    include("test_hertz.jl")
    include("test_traction.jl")
    include("test_nondim.jl")
    include("test_quaternion.jl")
    include("test_quasistatic.jl")
    include("test_dynamics.jl")
end
