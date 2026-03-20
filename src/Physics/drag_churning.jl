# =====================================================================
# Physics/drag_churning.jl — Lubricant drag & churning losses
#
# Port of dynamics_numba.py::_drag_force, _churning_moment, churning.py
# Generic types: supports Float64 and ForwardDiff.Dual
# =====================================================================

"""
    drag_coefficient_SN(Re) → C_D

Schiller-Naumann sphere drag: C_D = 24/Re + 6/(1+√Re) + 0.4.
Valid for Re < 2×10⁵.
"""
@inline function drag_coefficient_SN(Re)
    Re < 1e-10 && return zero(Re)
    return 24.0 / (Re + 1e-6) + 6.0 / (1.0 + sqrt(Re + 1e-6)) + 0.4
end

"""
    drag_force(V_abs, d, ρ, μ, web) → F_drag

Sphere drag force opposing orbital translation.
F = C_D · ½ρV² · A,  A = π/4·d² - w·d
"""
@inline function drag_force(V_abs, d, ρ, μ, web)
    T_ret = typeof(V_abs * ρ)
    (V_abs < 1e-15 || ρ < 1e-10) && return zero(T_ret)
    Re = ρ * V_abs * d / μ
    C_D = drag_coefficient_SN(Re)
    A = max(0.25 * π * d^2 - web * d, zero(T_ret))
    return C_D * 0.5 * ρ * V_abs^2 * A
end

"""
    churning_moment(ω_abs, r, ρ, μ) → M_churn

Disk-approximation churning moment (Eq 5.11).
M_e = ½ρω²r⁵·C_n, with laminar/turbulent C_n.
Generic: supports Float64 and ForwardDiff.Dual.
"""
function churning_moment(ω_abs, r, ρ, μ)
    T_ret = typeof(ω_abs * r * ρ)
    (ω_abs < 1e-10 || ρ < 1e-10 || r < 1e-10) && return zero(T_ret)
    Re = ρ * r^2 * ω_abs / μ
    Re < 1e-10 && return zero(T_ret)

    # Smooth laminar/turbulent transition (C∞) via blending, avoids if/else branch for AD
    C_lam = 3.87 / sqrt(Re + 1e-30)
    C_turb = 0.146 / (Re + 1e-30)^0.20
    σ = 1.0 / (1.0 + exp(-(Re - 300_000.0) / 30_000.0))  # sigmoid blend
    C_n = (1.0 - σ) * C_lam + σ * C_turb

    return 0.5 * ρ * ω_abs^2 * r^5 * C_n
end

"""
    effective_density(ρ_oil, ρ_air, fill) → ρ_eff
"""
@inline function effective_density(ρ_oil, ρ_air, fill)
    f = clamp(fill, 0.0, 1.0)
    return ρ_oil * f + ρ_air * (1.0 - f)
end
