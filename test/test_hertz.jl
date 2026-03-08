@testset "Hertz Contact" begin
    # 7210B bearing reference
    geom, mat = bearing_7210B()

    @testset "BearingGeometry construction" begin
        @test geom.d > 0
        @test geom.n_balls == 16
        @test r_ball(geom) ≈ 6.35e-3
        @test ball_mass(geom) > 0
        @test ball_inertia(geom) > 0
        @test ball_spacing(geom) ≈ 2π / 16
    end

    @testset "MaterialParams" begin
        @test mat.E ≈ 2.08e11
        @test mat.nu ≈ 0.3
        @test E_star(mat) > 0
        E_c = composite_modulus(mat, mat)
        @test E_c > 0
        # For same material: E* = E / (2(1-ν²))
        @test E_c ≈ E_star(mat)
    end

    @testset "Curvature sums" begin
        Σρ_i = sum_rho_inner(geom)
        Σρ_o = sum_rho_outer(geom)
        @test Σρ_i > 0   # inner must be positive (convex dominant)
        @test Σρ_o > 0   # outer also positive but smaller
        @test Σρ_i > Σρ_o # inner contact is tighter
    end

    @testset "Hertz stiffness computation" begin
        E_c = composite_modulus(mat, mat)
        h_inner, h_outer = create_bearing_hertz(geom, mat)

        @test h_inner.Upsilon > 0
        @test h_outer.Upsilon > 0
        @test h_inner.k > 1      # elliptic contact, a > b
        @test h_outer.k > 1

        # Contact load at 10μm penetration
        δ = 10e-6
        Q_inner = hertz_contact_load(δ, h_inner.Upsilon)
        Q_outer = hertz_contact_load(δ, h_outer.Upsilon)
        @test Q_inner > 0
        @test Q_outer > 0
        @test Q_inner > Q_outer  # inner is stiffer
    end

    @testset "Smooth Hertz transition" begin
        @test smooth_hertz_delta(1e-5) ≈ 1e-5 atol=1e-10
        @test smooth_hertz_delta(-1e-5) ≈ 0.0 atol=1e-10
        @test smooth_hertz_delta(0.0) > 0   # tiny, not zero
    end

    @testset "F(κ) monotonicity" begin
        # F should be monotonically increasing for k > 1
        ks = [1.1, 2.0, 5.0, 10.0, 50.0]
        Fs = [hertz_F_of_k(k) for k in ks]
        for i in 1:length(Fs)-1
            @test Fs[i] < Fs[i+1]
        end
        @test Fs[1] > 0
        @test Fs[end] < 1.0
    end
end
