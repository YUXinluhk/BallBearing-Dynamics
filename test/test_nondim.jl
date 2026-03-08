@testset "Non-dimensionalization" begin
    geom, mat = bearing_7210B()
    h_inner, h_outer = create_bearing_hertz(geom, mat)
    s = Scales(geom, h_inner)

    @testset "Scale magnitudes" begin
        @test s.Q > 0
        @test s.L ≈ r_ball(geom)
        @test s.T > 0
        @test s.V ≈ s.L / s.T
        @test s.W ≈ 1 / s.T
    end

    @testset "Roundtrip conversion" begin
        F = 1234.5  # [N]
        @test dim_force(s, nondim_force(s, F)) ≈ F atol=1e-10

        x = 0.001  # [m]
        @test dim_length(s, nondim_length(s, x)) ≈ x atol=1e-15

        t = 0.01  # [s]
        @test dim_time(s, nondim_time(s, t)) ≈ t atol=1e-15

        v = 5.0  # [m/s]
        @test dim_vel(s, nondim_vel(s, v)) ≈ v atol=1e-12

        ω = 1000.0  # [rad/s]
        @test dim_angvel(s, nondim_angvel(s, ω)) ≈ ω atol=1e-10
    end

    @testset "Inner Hertz Y★ ≈ 1" begin
        Y_star = nondim_stiffness(h_inner, s)
        # By construction of scales, Y_i★ should be ≈ 1
        @test Y_star ≈ 1.0 atol=0.01
    end
end
