# =====================================================================
# Transforms/frames.jl — Azimuth, contact, and reference frame rotations
#
# Port of transforms.py::T_azimuth, T_ac, contact angle helpers
# =====================================================================

using StaticArrays

"""
    T_azimuth(θ) → SMatrix{3,3}

Rotation from bearing coordinates to azimuth-plane coordinates.
X stays axial, Y/Z rotate by azimuth angle θ.
"""
@inline function T_azimuth(θ::Float64)
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
Azimuth-plane basis ordering is `[x, θ, r]`.
α₁ = primary axial-radial contact angle.
α₂ = lateral tilt of the contact normal out of the x-r plane.
"""
@inline function T_ac(α₁::Float64, α₂::Float64)
    s1, c1 = sincos(α₁)
    s2, c2 = sincos(α₂)
    @SMatrix [
        s1 * c2   s2    c1 * c2;
        -s1 * s2  c2    -c1 * s2;
        -c1       0.0   s1
    ]
end

"""
    contact_angles_from_direction(dx, dθ, dr) → (α₁, α₂)

Given the 3D displacement from ball center to race groove center expressed in the
local azimuth basis `[x, θ, r]`, compute the primary contact angle `α₁` and the
lateral tilt angle `α₂`.
"""
@inline function contact_angles_from_direction(dx::Float64, dθ::Float64, dr::Float64)
    α₁ = atan(dx, dr)
    L = sqrt(dx^2 + dθ^2 + dr^2 + 1e-30)
    α₂ = asin(clamp(dθ / L, -0.9999, 0.9999))
    return (α₁, α₂)
end

@inline function contact_angles_from_direction(dx::Float64, dr::Float64)
    return contact_angles_from_direction(dx, 0.0, dr)
end

"""
    contact_basis_from_angles(α₁, α₂) → (n_ac, t_lat_ac, t_roll_ac)

Return the orthonormal contact basis vectors in the local azimuth basis
`[x, θ, r]`. These correspond to the rows of `T_ac(α₁, α₂)`.
"""
@inline function contact_basis_from_angles(α₁::Float64, α₂::Float64)
    Tac = T_ac(α₁, α₂)
    n_ac = SVector{3,Float64}(Tac[1, 1], Tac[1, 2], Tac[1, 3])
    t_lat_ac = SVector{3,Float64}(Tac[2, 1], Tac[2, 2], Tac[2, 3])
    t_roll_ac = SVector{3,Float64}(Tac[3, 1], Tac[3, 2], Tac[3, 3])
    return (n_ac, t_lat_ac, t_roll_ac)
end

"""
    race_groove_center_pos(geom, θ, is_inner) → (x_gc, r_gc)

Race groove curvature center in azimuth-plane coordinates.
For initial geometry: x = 0, r = D_i/2 or D_o/2.
"""
function race_groove_center_pos(geom::BearingGeometry, θ::Float64, is_inner::Bool)
    if is_inner
        r_gc = D_i(geom) / 2.0
    else
        r_gc = D_o(geom) / 2.0
    end
    x_gc = 0.0  # groove center on midplane by default
    return (x_gc, r_gc)
end

"""
    inertial_from_cylindrical(x, r, θ) → SVector{3}

Convert cylindrical (x, r, θ) to Cartesian inertial (X, Y, Z).
X = x, Y = -r·sin(θ), Z = r·cos(θ).
"""
@inline function inertial_from_cylindrical(x::Float64, r::Float64, θ::Float64)
    sθ, cθ = sincos(θ)
    SVector{3,Float64}(x, -r * sθ, r * cθ)
end

"""
    velocity_inertial_from_cylindrical(ẋ, ṙ, θ̇, r, θ) → SVector{3}

Cartesian velocity from cylindrical time derivatives.
"""
@inline function velocity_inertial_from_cylindrical(ẋ, ṙ, θ̇, r, θ)
    sθ, cθ = sincos(θ)
    SVector{3,Float64}(ẋ, -ṙ*sθ - r*θ̇*cθ, ṙ*cθ - r*θ̇*sθ)
end
