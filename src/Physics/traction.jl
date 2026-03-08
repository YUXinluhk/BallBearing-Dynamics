# =====================================================================
# Physics/traction.jl — 4-parameter traction curve + spin moment
#
# Port of dynamics_numba.py::_traction_coeff, _spin_moment
# Generic types: supports Float64 and ForwardDiff.Dual
# =====================================================================

"""
    traction_coefficient(u_abs, A, B, C, D) → κ

4-parameter SKF traction curve: κ = (A + B·u)·exp(-C·u) + D.
Includes micro-gradient (1e-6·u) to stabilise Jacobian.
"""
@inline function traction_coefficient(u_abs, A, B, C, D)
    κ = (A + B * u_abs) * exp(-C * u_abs) + D + 1e-6 * u_abs
    return κ
end

"""
    tehd_traction_coefficient(u_abs, P_mean, a_hz, u_entrain,
                              A, B, C, D, Lambda_LSS, beta_temp) → κ_tehd

Lumped thermo-elasto-hydrodynamic (TEHD) traction model embedding:
1. Macro isothermal 4-parameter traction (SKF curve)
2. Bair-Winer LSS truncation (branchless `hypot`, C∞)
3. Archard contact flash temperature (ξ_steel ≈ 12000)
4. Nahme viscosity thermal softening
5. Elastic wave speed limiting (C_shear ≈ 3200 m/s)

This closes the Harris slip-heat paradox: at high slide velocities,
flash temperature → viscosity drop → traction self-caps → bounded Jacobian.
"""
@inline function tehd_traction_coefficient(u_abs, P_mean, a_hz, u_entrain,
    A, B, C, D, Lambda_LSS, beta_temp)
    eps_v = 1e-6
    u_s_smooth = abs(u_abs) + eps_v

    # 1. Macro isothermal traction coefficient (4-parameter model)
    mu_iso = (A + B * u_s_smooth) * exp(-C * u_s_smooth) + D

    # Short-circuit if no contact pressure
    (P_mean <= eps_v) && return mu_iso

    # 2. Bair-Winer LSS truncation: τ_max = Λ_LSS · P_mean
    #    μ_eff = μ_iso / √(1 + (μ_iso/Λ)²)  — branchless, C∞ smooth
    mu_eff = mu_iso / hypot(1.0, mu_iso / (Lambda_LSS + eps_v))

    # 3. Archard contact flash temperature
    #    q = μ·P·Δu,  ΔT_flash = 1.11·q·√(a/2) / (ξ·√|u_e|)
    #    ξ_steel ≈ 12000 W·s^0.5/(m²·K) (thermal inertia)
    q_heat = mu_eff * P_mean * u_s_smooth
    dT_flash = (1.11 * q_heat * sqrt(0.5 * abs(a_hz) + eps_v)) /
               (12000.0 * sqrt(abs(u_entrain) + eps_v))

    # 4. Nahme thermal softening: μ → μ / √(1 + β·ΔT)
    phi_T = 1.0 / sqrt(1.0 + beta_temp * dT_flash)

    # 5. Elastic wave speed limit (C_shear ≈ 3200 m/s for steel)
    psi_wave = exp(-(u_s_smooth / 3200.0)^2)

    return mu_eff * phi_T * psi_wave
end

"""
    spin_moment(ω_spin, Q, a, b, μ_spin, E2) → M_spin

Harris spin friction moment with precomputed elliptic integral.
M_spin = (3/8)·μ·Q·a·E2 · tanh(ω/ε)
E2 = complete elliptic integral E(m) precomputed from contact geometry.
"""
function spin_moment(ω_spin, Q, a, b, μ_spin, E2)
    Q < 1e-15 && return zero(Q)
    a < 1e-15 && return zero(a)

    # Harris formula: M_s = (3/8) μ Q a E(m)
    M = (3.0 / 8.0) * μ_spin * Q * a * E2

    # C∞ sign function: tanh(ω_spin/ε_spin)
    ε_spin = 1.0  # [rad/s] transition width
    sign_smooth = tanh(ω_spin / ε_spin)

    return M * sign_smooth
end
