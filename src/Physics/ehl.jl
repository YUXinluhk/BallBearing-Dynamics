# =====================================================================
# Physics/ehl.jl — Elastohydrodynamic Lubrication film thickness
#
# Port of ehl_traction.py film thickness and viscosity functions.
# Stateless pure functions: (params_in) → (h, ϕ_T)
# Generic types: supports Float64 and ForwardDiff.Dual
# =====================================================================

"""
    viscosity(p, T, lub) → μ

Lubricant dynamic viscosity at given pressure and temperature.
μ = μ₀ · exp(α·p + β·(T₀ - T))
"""
@inline function viscosity(p, T, lub::LubricantParams)
    exponent = lub.alpha_pv * p + lub.beta_temp * (lub.T_0 - T)
    exponent = clamp(exponent, -30.0, 30.0)  # overflow protection
    return lub.mu_0 * exp(exponent)
end

"""
    thermal_correction(μ₀, u_mean, K_th, T₀) → ϕ_T

Thermal correction factor for film thickness.
ϕ_T ≈ 1/(1 + 0.1·Q_m^{0.64}), Q_m = 2μ₀u²/(K_th·T₀)
"""
@inline function thermal_correction(μ₀, u_mean, K_th, T₀)
    u² = u_mean^2
    T_ret = typeof(μ₀ * u²)
    (u² < 1e-30 || K_th < 1e-30 || T₀ < 1e-10) && return one(T_ret)
    Q_m = 2.0 * μ₀ * u² / (K_th * T₀)
    return 1.0 / (1.0 + 0.1 * Q_m^0.64)
end

"""
    film_thickness_hd(u_mean, Q, E_prime, R_eff, κ_e, μ₀, α_pv, K_th, T₀) → h

Hamrock-Dowson central film thickness with thermal correction.
h/R = 2.69·U^{0.67}·G^{0.53}·(1 - 0.61·exp(-0.73κ))/W^{0.067}
"""
function film_thickness_hd(u_mean, Q, E_prime, R_eff, κ_e,
    μ₀, α_pv, K_th, T₀)
    T_ret = typeof(u_mean * Q)
    (abs(u_mean) < 1e-15 || E_prime < 1e-10 || R_eff < 1e-15) && return zero(T_ret)

    E_2star = 2.0 * E_prime   # EHL convention E' = 2E*

    # Dimensionless groups
    U_param = μ₀ * abs(u_mean) / (E_2star * R_eff)
    G_param = α_pv * E_2star
    W_param = max(abs(Q) / (E_2star * R_eff^2), 1e-30)

    # Log-domain computation to avoid underflow
    ln_h_R = log(2.69) + 0.67 * log(U_param + 1e-30) +
             0.53 * log(G_param + 1e-30) - 0.067 * log(W_param)

    # Ellipticity correction
    ellip_corr = 1.0 - 0.61 * exp(-0.73 * κ_e)

    h = R_eff * exp(ln_h_R) * ellip_corr

    # Thermal correction
    ϕ_T = thermal_correction(μ₀, u_mean, K_th, T₀)
    h *= ϕ_T

    return max(h, zero(T_ret))
end
