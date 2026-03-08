# =====================================================================
# Physics/rolling_resistance.jl — Palmgren viscous rolling resistance
#
# Port of rolling_resistance.py
# =====================================================================

"""
    palmgren_Mrr(d_m_mm, f₀, μ_Pa_s, speed_rpm, ρ_lub) → M_rr [N·mm]

Palmgren rolling resistance torque.
M_rr = f₀ × 1e-7 × (ν·n)^{2/3} × d_m³
"""
function palmgren_Mrr(d_m_mm::Float64, f₀::Float64,
                      μ_Pa_s::Float64, speed_rpm::Float64,
                      ρ_lub::Float64=860.0)
    ν_mm² = μ_Pa_s / ρ_lub * 1e6   # → mm²/s
    νn = ν_mm² * speed_rpm

    if νn < 2000.0
        return 160e-7 * d_m_mm^3
    else
        return f₀ * 1e-7 * νn^(2/3) * d_m_mm^3
    end
end

"""
    palmgren_power(d_m_mm, f₀, μ_Pa_s, speed_rpm, ρ_lub) → P_rr [W]
"""
function palmgren_power(d_m_mm::Float64, f₀::Float64,
                        μ_Pa_s::Float64, speed_rpm::Float64,
                        ρ_lub::Float64=860.0)
    M_rr = palmgren_Mrr(d_m_mm, f₀, μ_Pa_s, speed_rpm, ρ_lub)
    ω = speed_rpm * π / 30.0
    return M_rr * 1e-3 * ω   # N·mm → N·m, × rad/s → W
end
