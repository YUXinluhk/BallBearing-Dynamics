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
α₁ = primary contact angle, α₂ = lateral tilt.
"""
@inline function T_ac(α₁::Float64, α₂::Float64)
    s1, c1 = sincos(α₁)
    s2, c2 = sincos(α₂)
    @SMatrix [
        s1*c2    s1*s2    c1;
        s2       -c2      0.0;
        -c1*c2   -c1*s2   s1
    ]
end

"""
    contact_angles_from_direction(dx, dr) → (α₁, α₂)

Given displacement from ball center to race groove center (axial dx, radial dr),
compute primary contact angle α₁.
"""
@inline function contact_angles_from_direction(dx::Float64, dr::Float64)
    α₁ = atan(dx, dr)   # atan2(axial, radial)
    return (α₁, 0.0)     # α₂ = 0 for axial symmetry
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
