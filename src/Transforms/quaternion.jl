# =====================================================================
# Transforms/quaternion.jl — Unit quaternion kinematics
#
# Replaces Euler-angle machinery (transforms.py::TrigCache, T_ib, C_inv).
# Uses Quaternions.jl for core operations.
# =====================================================================

using Quaternions: Quaternion, QuaternionF64, quat, normalize, imag_part
using StaticArrays

"""Unit quaternion identity"""
const QUAT_IDENTITY = QuaternionF64(1.0, 0.0, 0.0, 0.0)

"""
    quat_from_components(s, v1, v2, v3) → Quaternion

Create quaternion q = s + v1·i + v2·j + v3·k.
"""
@inline quat_from_components(s, v1, v2, v3) = QuaternionF64(s, v1, v2, v3)

"""
    quat_derivative(q, ω_body) → dq/dt

Quaternion kinematic equation: q̇ = ½ q ⊗ [0, ω_body].
"""
@inline function quat_derivative(q::QuaternionF64, ω::SVector{3,Float64})
    ω_q = QuaternionF64(0.0, ω[1], ω[2], ω[3])
    return 0.5 * q * ω_q
end

"""
    rotate_vector(q, v) → v'

Rotate vector v by quaternion q: v' = q v q*.
"""
@inline function rotate_vector(q::QuaternionF64, v::SVector{3,Float64})
    v_q = QuaternionF64(0.0, v[1], v[2], v[3])
    r = q * v_q * conj(q)
    return SVector{3,Float64}(imag_part(r)...)
end

"""
    inv_rotate_vector(q, v) → v'

Inverse-rotate vector v by quaternion q: v' = q* v q.
"""
@inline function inv_rotate_vector(q::QuaternionF64, v::SVector{3,Float64})
    v_q = QuaternionF64(0.0, v[1], v[2], v[3])
    r = conj(q) * v_q * q
    return SVector{3,Float64}(imag_part(r)...)
end

"""
    quat_to_rotmat(q) → SMatrix{3,3}

Convert quaternion to 3×3 rotation matrix (body → inertial).
"""
function quat_to_rotmat(q::QuaternionF64)
    s = real(q)
    v1, v2, v3 = imag_part(q)

    # Normalize
    n = s^2 + v1^2 + v2^2 + v3^2
    n < 1e-30 && return SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    fac = 2.0 / n

    xx = v1 * v1 * fac; yy = v2 * v2 * fac; zz = v3 * v3 * fac
    xy = v1 * v2 * fac; xz = v1 * v3 * fac; yz = v2 * v3 * fac
    wx = s * v1 * fac;  wy = s * v2 * fac;  wz = s * v3 * fac

    @SMatrix [
        1.0 - yy - zz    xy - wz          xz + wy;
        xy + wz           1.0 - xx - zz    yz - wx;
        xz - wy           yz + wx          1.0 - xx - yy
    ]
end

"""
    quat_renormalize(q) → q_normalized

Enforce unit norm constraint: q → q / |q|.
"""
@inline function quat_renormalize(q::QuaternionF64)
    n = sqrt(real(q)^2 + sum(imag_part(q).^2))
    n < 1e-30 && return QUAT_IDENTITY
    return q / n
end

"""
    quat_from_euler_zyx(η, ξ, λ) → Quaternion

Convert Tait-Bryan Z-Y-X angles to quaternion.
Useful for converting initial conditions from the Python code.
"""
function quat_from_euler_zyx(η::Float64, ξ::Float64, λ::Float64)
    c1, s1 = cos(η/2), sin(η/2)
    c2, s2 = cos(ξ/2), sin(ξ/2)
    c3, s3 = cos(λ/2), sin(λ/2)

    w = c1*c2*c3 + s1*s2*s3
    x = s1*c2*c3 - c1*s2*s3
    y = c1*s2*c3 + s1*c2*s3
    z = c1*c2*s3 - s1*s2*c3

    return QuaternionF64(w, x, y, z)
end

"""
    omega_body_from_euler(η, ξ, λ, η̇, ξ̇, λ̇) → SVector{3}

Angular velocity in body frame from Euler angle rates.
(For converting initial conditions.)
"""
function omega_body_from_euler(η, ξ, λ, η̇, ξ̇, λ̇)
    cξ, sξ = cos(ξ), sin(ξ)
    cλ, sλ = cos(λ), sin(λ)
    ω_x = η̇ * cξ * cλ + ξ̇ * sλ
    ω_y = -η̇ * cξ * sλ + ξ̇ * cλ
    ω_z = η̇ * sξ + λ̇
    return SVector{3,Float64}(ω_x, ω_y, ω_z)
end
