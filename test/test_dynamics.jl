@testset "Dynamic Simulation" begin
    geom, mat = bearing_7210B()
    lub = lubricant_mil_l_23699()
    trac = traction_params_default()
    cage = cage_from_bearing(geom)

    @testset "State vector layout" begin
        Z = geom.n_balls
        @test n_state(Z) == 18 + 13Z
        @test n_state(16) == 226

        u = zeros(n_state(Z))
        # Check views don't overlap
        @test length(ir_pos_view(u, Z)) == N_IR_POS
        @test length(ball_pos_view(u, 1, Z)) == N_BALL_POS
        @test length(cage_vel_view(u, Z)) == N_CAGE_VEL
    end

    @testset "Quaternion in state vector" begin
        Z = geom.n_balls
        u = zeros(n_state(Z))
        # Set identity quaternion for ball 1
        set_ball_quat!(u, 1, Z, QUAT_IDENTITY)
        q = ball_quat(u, 1, Z)
        @test real(q) ≈ 1.0
        @test all(imag_part(q) .≈ 0.0)
    end

    @testset "Parameter vector construction" begin
        h_inner, h_outer = create_bearing_hertz(geom, mat)
        s = Scales(geom, h_inner)
        config = SimulationConfig(
            t_end=0.001, dt_output=1e-5,
            inner_race_speed=1000.0, F_axial=2000.0
        )
        params = build_params(geom, mat, lub, trac, cage, config, s, h_inner, h_outer)

        @test length(params) == N_PARAMS
        @test params[P_NBALL] ≈ 16.0
        @test params[P_YI] > 0
        @test params[P_FA] > 0
    end
end
