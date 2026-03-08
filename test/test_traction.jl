@testset "Traction Curve" begin
    tp = traction_params_default()

    @testset "Boundary conditions" begin
        # κ(0) ≈ kappa_0
        κ_0 = traction_coefficient(0.0, tp.A, tp.B, tp.C, tp.D)
        @test κ_0 ≈ tp.kappa_0 atol=1e-3

        # κ(∞) → kappa_inf
        κ_inf = traction_coefficient(100.0, tp.A, tp.B, tp.C, tp.D)
        @test κ_inf ≈ tp.kappa_inf atol=5e-3
    end

    @testset "Peak location" begin
        # Find peak numerically
        us = range(0.0, 5.0, length=1000)
        κs = [traction_coefficient(u, tp.A, tp.B, tp.C, tp.D) for u in us]
        u_peak = us[argmax(κs)]
        @test u_peak ≈ tp.u_m atol=0.1  # within 20%
    end

    @testset "Spin moment" begin
        # No load → no spin moment
        M0 = spin_moment(100.0, 0.0, 1e-3, 0.5e-3, 0.06)
        @test M0 ≈ 0.0

        # With load, spin moment opposes spin
        M_pos = spin_moment(100.0, 1000.0, 1e-3, 0.5e-3, 0.06)
        M_neg = spin_moment(-100.0, 1000.0, 1e-3, 0.5e-3, 0.06)
        @test M_pos > 0
        @test M_neg < 0
        @test abs(M_pos) ≈ abs(M_neg) atol=1e-6
    end
end
