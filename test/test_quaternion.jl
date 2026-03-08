using StaticArrays

@testset "Quaternion Kinematics" begin

    @testset "Identity" begin
        q = QUAT_IDENTITY
        v = SVector{3,Float64}(1.0, 2.0, 3.0)
        v_rot = rotate_vector(q, v)
        @test v_rot ≈ v atol=1e-14
    end

    @testset "90° rotation about Z" begin
        # q = cos(45°) + sin(45°)·k
        q = QuaternionF64(cos(π/4), 0.0, 0.0, sin(π/4))
        v = SVector{3,Float64}(1.0, 0.0, 0.0)
        v_rot = rotate_vector(q, v)
        @test v_rot[1] ≈ 0.0 atol=1e-14
        @test v_rot[2] ≈ 1.0 atol=1e-14
        @test v_rot[3] ≈ 0.0 atol=1e-14
    end

    @testset "Inverse rotation" begin
        q = QuaternionF64(cos(π/6), sin(π/6), 0.0, 0.0)
        v = SVector{3,Float64}(1.0, 1.0, 1.0)
        v_rot = rotate_vector(q, v)
        v_back = inv_rotate_vector(q, v_rot)
        @test v_back ≈ v atol=1e-13
    end

    @testset "Quaternion derivative consistency" begin
        q = quat_renormalize(QuaternionF64(1.0, 0.1, 0.2, 0.3))
        ω = SVector{3,Float64}(10.0, 0.0, 0.0)
        dq = quat_derivative(q, ω)

        # dq should be orthogonal to q: Re(q* · dq) ≈ 0
        # (since d|q|²/dt = 2Re(q*·dq) = 0 for unit quaternions)
        orth = real(conj(q) * dq)
        @test abs(orth) < 1e-10
    end

    @testset "Renormalization" begin
        q = QuaternionF64(10.0, 0.1, 0.2, 0.3)
        q_n = quat_renormalize(q)
        n = sqrt(sum(imag_part(q_n).^2) + real(q_n)^2)
        @test n ≈ 1.0 atol=1e-15
    end

    @testset "Euler → quaternion roundtrip" begin
        η, ξ, λ = 0.1, 0.2, 0.3
        q = quat_from_euler_zyx(η, ξ, λ)
        R = quat_to_rotmat(q)
        # Check that R is orthogonal
        @test norm(R * R' - I) < 1e-12
    end
end
