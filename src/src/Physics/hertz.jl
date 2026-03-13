<<<<<<< Current (Your changes)
=======
# =====================================================================
# Physics/hertz.jl — Hertz point-contact elasticity solution
#
# Port of hertz.py. Uses Elliptic.jl for exact K(m), E(m).
# Single source of truth for both quasi-static and dynamics.
# All functions are STATELESS pure functions.
# =====================================================================

using Elliptic: K as ellipK, E as ellipE

# ── Core: F(k) and bisection ─────────────────────────────────────────

"""
    hertz_F_of_k(k) → F

F(k) = [(k²+1)E(m) - 2K(m)] / [(k²-1)E(m)],  m = 1 - 1/k².
Monotonically increasing for k > 1; F(1)=0, F(∞)→1.
"""
@inline function hertz_F_of_k(k::Float64)
    k <= 1.0 && return 0.0
    k² = k * k
    m = 1.0 - 1.0 / k²
    K_val = ellipK(m)
    E_val = ellipE(m)
    denom = (k² - 1.0) * E_val
    abs(denom) < 1e-30 && return 0.0
    return ((k² + 1.0) * E_val - 2.0 * K_val) / denom
end

"""
    solve_kappa_bisection(F_rho) → (κ, K_val, E_val)

Solve F(κ) = F_rho via bisection with adaptive upper bound.
"""
function solve_kappa_bisection(F_rho::Float64)
    F_rho < 1e-10 && return (1.0, Float64(π) / 2, Float64(π) / 2)
    F_rho_c = min(F_rho, 0.9999)

    k_lo = 1.0 + 1e-8
    k_hi = 2.0
    # Expand upper bound until it brackets solution
    for _ in 1:50
        hertz_F_of_k(k_hi) >= F_rho_c && break
        k_hi *= 2.0
    end
    # Bisection (40 iterations → ~1e-12 precision)
    for _ in 1:40
        k_mid = 0.5 * (k_lo + k_hi)
        if hertz_F_of_k(k_mid) < F_rho_c
            k_lo = k_mid
        else
            k_hi = k_mid
        end
    end
    κ = 0.5 * (k_lo + k_hi)
    m = 1.0 - 1.0 / (κ * κ)
    return (κ, ellipK(m), ellipE(m))
end

# ── Canonical Hertz computation ──────────────────────────────────────

"""
    hertz_Y_and_ab(cos_alpha, d, d_m, f, E_prime, is_inner)
        → (Y, a_star, b_star, sum_rho)

Canonical Hertz stiffness computation. Y gives Q = Y·δ^{3/2}.

Exactly matches Python `_hertz_Y_and_ab` (hertz.py L95-144).
"""
function hertz_Y_and_ab(cos_alpha::Float64, d::Float64, d_m::Float64,
    f::Float64, E_prime::Float64, is_inner::Bool)
    ca = clamp(cos_alpha, 1e-4, 0.9999)
    ρ_ball = 2.0 / d

    if is_inner
        ρ_race_roll = 2.0 * ca / (d_m - d * ca)
    else
        ρ_race_roll = -2.0 * ca / (d_m + d * ca)
    end
    ρ_race_cross = -1.0 / (f * d)
    Σρ = 2.0 * ρ_ball + ρ_race_roll + ρ_race_cross

    Σρ < 1e-20 && return (0.0, 1.0, 1.0, Σρ)

    F_ρ = min(abs(ρ_race_cross - ρ_race_roll) / Σρ, 0.9999)

    κ, K_val, E_val = solve_kappa_bisection(F_ρ)

    factor = π / (2.0 * κ^2 * E_val)
    factor < 1e-30 && return (0.0, 1.0, 1.0, Σρ)

    a_star = (2.0 * κ^2 * E_val / π)^(1 / 3)
    b_star = (2.0 * E_val / (π * κ))^(1 / 3)
    δ_star = (2.0 * K_val / π) * factor^(1 / 3)

    δ_star < 1e-30 && return (0.0, a_star, b_star, Σρ)

    # Key formula: Y = 4√2·E′ / (3√Σρ · δ★^{3/2})
    Y = 4.0 * √2.0 * E_prime / (3.0 * √Σρ * δ_star^1.5)

    return (Y, a_star, b_star, Σρ)
end

"""
    hertz_runtime_contact(cos_alpha, d, d_m, f, E_prime, is_inner)
        → (Y, a_star, b_star, sum_rho, E2)

Runtime Hertz contact constants for the instantaneous contact angle.
Returns the complete elliptic integral `E2 = E(m)` as well, so dynamics can
update the Harris spin moment consistently with the current contact ellipse.
"""
function hertz_runtime_contact(cos_alpha::Float64, d::Float64, d_m::Float64,
    f::Float64, E_prime::Float64, is_inner::Bool)
    ca = clamp(cos_alpha, 1e-4, 0.9999)
    ρ_ball = 2.0 / d

    if is_inner
        ρ_race_roll = 2.0 * ca / (d_m - d * ca)
    else
        ρ_race_roll = -2.0 * ca / (d_m + d * ca)
    end
    ρ_race_cross = -1.0 / (f * d)
    Σρ = 2.0 * ρ_ball + ρ_race_roll + ρ_race_cross

    Σρ < 1e-20 && return (0.0, 1.0, 1.0, Σρ, Float64(π) / 2)

    F_ρ = min(abs(ρ_race_cross - ρ_race_roll) / Σρ, 0.9999)
    κ, K_val, E_val = solve_kappa_bisection(F_ρ)

    factor = π / (2.0 * κ^2 * E_val)
    factor < 1e-30 && return (0.0, 1.0, 1.0, Σρ, E_val)

    a_star = (2.0 * κ^2 * E_val / π)^(1 / 3)
    b_star = (2.0 * E_val / (π * κ))^(1 / 3)
    δ_star = (2.0 * K_val / π) * factor^(1 / 3)

    δ_star < 1e-30 && return (0.0, a_star, b_star, Σρ, E_val)

    Y = 4.0 * √2.0 * E_prime / (3.0 * √Σρ * δ_star^1.5)
    return (Y, a_star, b_star, Σρ, E_val)
end

"""
    hertz_from_curvatures(sum_rho, F_rho, E_prime)
        → (κ, K_val, E_val, a★, b★, δ★, Y)

Hertz solution from curvature parameters (used by quasi-static solver).
"""
function hertz_from_curvatures(Σρ::Float64, F_ρ::Float64, E_prime::Float64)
    κ, K_val, E_val = solve_kappa_bisection(F_ρ)

    a_star = (2.0 * κ^2 * E_val / π)^(1 / 3)
    b_star = (2.0 * E_val / (π * κ))^(1 / 3)
    δ_star = (2.0 * K_val / π) * (π / (2.0 * κ^2 * E_val))^(1 / 3)

    Y = if Σρ > 0.0 && δ_star > 0.0
        4.0 * √2.0 * E_prime / (3.0 * √Σρ * δ_star^1.5)
    else
        0.0
    end

    return (κ, K_val, E_val, a_star, b_star, δ_star, Y)
end

"""
    hertz_from_curvatures_vec(sum_rho_arr, F_rho_arr, E_prime)

Vectorized Hertz solver returning NamedTuple of arrays.
"""
function hertz_from_curvatures_vec(Σρ_arr::Vector{Float64},
    F_ρ_arr::Vector{Float64},
    E_prime::Float64)
    n = length(Σρ_arr)
    κ = Vector{Float64}(undef, n)
    K1 = Vector{Float64}(undef, n)
    E2 = Vector{Float64}(undef, n)
    a_s = Vector{Float64}(undef, n)
    b_s = Vector{Float64}(undef, n)
    δ_s = Vector{Float64}(undef, n)
    Y = Vector{Float64}(undef, n)

    for j in 1:n
        κ[j], K1[j], E2[j], a_s[j], b_s[j], δ_s[j], Y[j] =
            hertz_from_curvatures(Σρ_arr[j], F_ρ_arr[j], E_prime)
    end
    return (kappa=κ, K1=K1, E2=E2, a_star=a_s, b_star=b_s, delta_star=δ_s, Upsilon=Y)
end

# ── Runtime contact functions ────────────────────────────────────────

"""
    hertz_contact_load(δ, Y) → Q

Q = Y · max(0, δ)^{3/2}
"""
@inline function hertz_contact_load(δ::Float64, Y::Float64)
    δ <= 0.0 && return 0.0
    return Y * δ^1.5
end

"""
    hertz_ab(Q, a★, b★, E′, Σρ) → (a, b)

Dimensional contact ellipse semi-axes from load.
"""
@inline function hertz_ab(Q::Float64, a_star::Float64, b_star::Float64,
    E_prime::Float64, Σρ::Float64)
    Q <= 0.0 && return (0.0, 0.0)
    c = (3Q / (2E_prime * Σρ))^(1 / 3)
    return (a_star * c, b_star * c)
end

"""
    smooth_hertz_delta(δ; ε=1e-12)

C∞-smooth transition: δ_sm = ½(δ + √(δ² + ε²)) ≈ max(0, δ).
"""
@inline function smooth_hertz_delta(δ; ε=1e-5)  # 【Bug7修复】ε^2=1e-10 survives Float64 at δ~1e-4
    return 0.5 * (δ + sqrt(δ^2 + ε^2))
end

# ── Precomputed contact struct ───────────────────────────────────────

"""
    HertzContact — immutable precomputed contact constants for a ball/race pair.
"""
struct HertzContact
    k::Float64          # contact ellipse ratio a/b
    K1::Float64         # complete elliptic integral K(m)
    E2::Float64         # complete elliptic integral E(m)
    a_star::Float64     # dimensionless semi-major axis
    b_star::Float64     # dimensionless semi-minor axis
    delta_star::Float64 # dimensionless deformation
    Upsilon::Float64    # stiffness: Q = Υ·δ^{3/2} [N/m^1.5]
    sum_rho::Float64    # Σρ̄ [1/m]
    E_prime::Float64    # composite modulus E* [Pa]
end

function HertzContact(; sum_rho::Float64, F_rho::Float64, E_prime::Float64)
    κ, K1, E2, a_s, b_s, δ_s, Y = hertz_from_curvatures(sum_rho, F_rho, E_prime)
    HertzContact(κ, K1, E2, a_s, b_s, δ_s, Y, sum_rho, E_prime)
end

"Create inner and outer HertzContact from BearingGeometry + MaterialParams"
function create_bearing_hertz(geom::BearingGeometry, mat::MaterialParams)
    E_p = composite_modulus(mat, mat)
    h_inner = HertzContact(sum_rho=sum_rho_inner(geom), F_rho=F_rho_inner(geom), E_prime=E_p)
    h_outer = HertzContact(sum_rho=sum_rho_outer(geom), F_rho=F_rho_outer(geom), E_prime=E_p)
    return (h_inner, h_outer)
end

"Non-dimensional stiffness: Υ★ = Υ·L^{3/2}/Q₀"
nondim_stiffness(h::HertzContact, s) = h.Upsilon * s.L^1.5 / s.Q
>>>>>>> Incoming (Background Agent changes)
