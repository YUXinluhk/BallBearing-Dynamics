# =====================================================================
# Transforms/frames.jl — Azimuth, contact, and reference frame rotations
#
# Port of transforms.py::T_azimuth, T_ac, contact angle helpers
# =====================================================================

using StaticArrays

"""
    T_azimuth(θ) → SMatrix{3,3}

Rotation from bearing coordinates to azimuth-plane coordinates.
Basis ordering: `[x, r, θ]` (matches kernel inline convention).
X stays axial, Y/Z rotate by azimuth angle θ.
"""
@inline function T_azimuth(θ)
    sθ, cθ = sincos(θ)
    @SMatrix [
        1.0   0.0   0.0;
        0.0  -sθ    cθ;
        0.0  -cθ   -sθ
    ]
end

"""
    T_ac(α₁, α₂) → SMatrix{3,3}

Rotation from azimuth plane to contact plane.
Azimuth-plane basis ordering is `[x, r, θ]`.
α₁ = primary axial-radial contact angle (atan(dx, dr)).
α₂ = lateral tilt of the contact normal out of the x-r plane.

Rows: [n_contact; t_lateral; t_rolling]
"""
@inline function T_ac(α₁, α₂)
    s1, c1 = sincos(α₁)
    s2, c2 = sincos(α₂)
    z = zero(typeof(s1))
    @SMatrix [
        s1 * c2   c1 * c2   s2;
       -c1        s1        z;
       -s1 * s2  -c1 * s2   c2
    ]
end

"""
    contact_basis_from_displacements(dx, dr, dθ)

【极速优化版】: 纯代数推导接触系基向量，彻底消灭 atan, asin, sin, cos！
相比基于角度的旋转矩阵，速度提升极高，且完全免疫 AD 奇异点。
"""
@inline function contact_basis_from_displacements(dx, dr, dθ)
    L_plane = sqrt(dx^2 + dr^2 + 1e-30)
    L_3d = sqrt(dx^2 + dr^2 + dθ^2 + 1e-30)
    
    s1 = dx / L_plane
    c1 = dr / L_plane
    s2 = dθ / L_3d
    c2 = L_plane / L_3d  # 纯代数推导，避免 clamp 和 sqrt
    
    T = typeof(s1 * s2)
    z = zero(T)
    
    n_ac      = SVector{3,T}(s1 * c2,   c1 * c2,   s2)
    t_lat_ac  = SVector{3,T}(-c1,       s1,        z)
    t_roll_ac = SVector{3,T}(-s1 * s2, -c1 * s2,   c2)
    
    return n_ac, t_lat_ac, t_roll_ac
end

"""
    smooth_signed_saturate(x, limit)

Smooth odd saturation used to regularize algebraic closures without introducing
hard kinks into the Jacobian.
"""
@inline function smooth_signed_saturate(x, limit)
    return limit * tanh(x / limit)
end

"""
    conservative_contact_pass_time(a_hz, b_hz, v_pass)

Estimate contact residence time using a conservative equivalent ellipse length.
`hypot(a, b)` is smoother than `max(a, b)` and gives a longer traversal scale.
"""
@inline function conservative_contact_pass_time(a_hz, b_hz, v_pass)
    l_eff = sqrt(a_hz^2 + b_hz^2 + 1e-24)
    return l_eff / sqrt(v_pass^2 + 1e-12)
end

"""
    race_groove_center_pos(geom, θ, is_inner) → (x_gc, r_gc)

Race groove curvature center in azimuth-plane coordinates.
For initial geometry: x = 0, r = D_i/2 or D_o/2.
"""
function race_groove_center_pos(geom::BearingGeometry, θ, is_inner::Bool)
    if is_inner
        r_gc = D_i(geom) / 2.0
    else
        r_gc = D_o(geom) / 2.0
    end
    x_gc = 0.0
    return (x_gc, r_gc)
end

"""
    inertial_from_cylindrical(x, r, θ) → SVector{3}

Convert cylindrical (x, r, θ) to Cartesian inertial (X, Y, Z).
X = x, Y = -r·sin(θ), Z = r·cos(θ).
"""
@inline function inertial_from_cylindrical(x, r, θ)
    sθ, cθ = sincos(θ)
    T = typeof(x * sθ)
    SVector{3, T}(x, -r * sθ, r * cθ)
end

"""
    velocity_inertial_from_cylindrical(ẋ, ṙ, θ̇, r, θ) → SVector{3}

Cartesian velocity from cylindrical time derivatives.
"""
@inline function velocity_inertial_from_cylindrical(ẋ, ṙ, θ̇, r, θ)
    sθ, cθ = sincos(θ)
    T = typeof(ẋ * sθ)
    SVector{3, T}(ẋ, -ṙ*sθ - r*θ̇*cθ, ṙ*cθ - r*θ̇*sθ)
end
