@testset "Quasi-Static Solver" begin
    geom, mat = bearing_7210B()

    @testset "7210B converges under axial load" begin
        result = solve_quasi_static(geom, mat;
                                    F_a=2000.0, n_rpm=0.0)
        @test result.converged
        @test all(result.Q_inner .> 0)
        @test all(result.alpha_inner .> 0)
    end

    @testset "Contact angles increase with load" begin
        r1 = solve_quasi_static(geom, mat; F_a=1000.0)
        r2 = solve_quasi_static(geom, mat; F_a=5000.0)
        # Higher axial load → larger inner contact angle
        @test mean(r2.alpha_inner) > mean(r1.alpha_inner)
    end

    @testset "Static: all balls equally loaded" begin
        result = solve_quasi_static(geom, mat; F_a=2000.0, n_rpm=0.0)
        Q_mean = mean(result.Q_inner)
        # For pure axial, all balls should carry nearly equal load
        for j in 1:geom.n_balls
            @test abs(result.Q_inner[j] - Q_mean) / Q_mean < 0.01
        end
    end

    @testset "Combined load: radial variation" begin
        result = solve_quasi_static(geom, mat;
                                    F_a=2000.0, F_rz=1000.0, n_rpm=0.0)
        @test result.converged
        # Radial load creates non-uniform distribution
        @test maximum(result.Q_inner) > minimum(result.Q_inner)
    end
end
