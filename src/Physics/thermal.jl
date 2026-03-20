# =====================================================================
# Physics/thermal.jl — LPTN 4-node thermal network
#
# Port of thermal.py. All pure functions.
# =====================================================================

"""
    setup_thermal_capacitances!(tp, geom) → ThermalParams

Auto-estimate heat capacitances from geometry (if NaN).
Returns a NEW ThermalParams with capacitances filled.
"""
function setup_thermal_capacitances(tp::ThermalParams, geom::BearingGeometry)
    c = tp.c_steel
    C_ir = isnan(tp.C_ir) ? inner_race_mass(geom) * c : tp.C_ir
    C_or = isnan(tp.C_or) ? inner_race_mass(geom) * 1.5 * c : tp.C_or
    C_ball = isnan(tp.C_ball) ? geom.n_balls * ball_mass(geom) * c : tp.C_ball
    C_oil = isnan(tp.C_oil) ? 5e-6 * 860.0 * 2000.0 : tp.C_oil

    return ThermalParams(
        T_ambient = tp.T_ambient, T_init = tp.T_init,
        C_ir = C_ir, C_or = C_or, C_ball = C_ball, C_oil = C_oil,
        G_ir_ball = tp.G_ir_ball, G_or_ball = tp.G_or_ball,
        G_ball_oil = tp.G_ball_oil, G_or_amb = tp.G_or_amb, G_oil_amb = tp.G_oil_amb,
        CTE = tp.CTE, T_ref = tp.T_ref, eta_ball = tp.eta_ball, c_steel = tp.c_steel,
    )
end

"""
    lptn_step(state, Q_ir, Q_or, Q_ball, Q_oil, dt, tp) → LPTNState

Forward Euler thermal ODE step: Cᵢ·dTᵢ/dt = Qᵢ,gen - Σ Gᵢⱼ·(Tᵢ - Tⱼ)
"""
function lptn_step(state::LPTNState, Q_ir::Float64, Q_or::Float64,
                   Q_ball::Float64, Q_oil::Float64,
                   dt::Float64, tp::ThermalParams)
    T = as_array(state)
    T_amb = tp.T_ambient

    dT_ir = (Q_ir - tp.G_ir_ball * (T[1] - T[3])) / tp.C_ir
    dT_or = (Q_or - tp.G_or_ball * (T[2] - T[3]) - tp.G_or_amb * (T[2] - T_amb)) / tp.C_or
    dT_ball = (Q_ball + tp.G_ir_ball * (T[1] - T[3]) + tp.G_or_ball * (T[2] - T[3]) -
               tp.G_ball_oil * (T[3] - T[4])) / tp.C_ball
    dT_oil = (Q_oil + tp.G_ball_oil * (T[3] - T[4]) - tp.G_oil_amb * (T[4] - T_amb)) / tp.C_oil

    return LPTNState(
        T[1] + dt * dT_ir,
        T[2] + dt * dT_or,
        T[3] + dt * dT_ball,
        T[4] + dt * dT_oil,
    )
end

"""
    μ₀_at_temp(T_oil, lub) → μ₀

Oil viscosity at temperature: μ₀(T) = μ₀,ref · exp(β·(T₀ - T))
"""
@inline function μ₀_at_temp(T_oil::Float64, lub::LubricantParams)
    exponent = clamp(lub.beta_temp * (lub.T_0 - T_oil), -30.0, 30.0)
    return lub.mu_0 * exp(exponent)
end

"""
    δr_thermal(T_ir, tp, R_raceway) → Δr

Inner race thermal radial expansion: Δr = CTE·(T_ir - T_ref)·R
"""
@inline function δr_thermal(T_ir::Float64, tp::ThermalParams, R_raceway::Float64)
    return tp.CTE * (T_ir - tp.T_ref) * R_raceway
end

"""
    distribute_heat(heat_avg, tp) → (Q_ir, Q_or, Q_ball, Q_oil)

Distribute average heat generation to LPTN nodes.
"""
function distribute_heat(heat_avg::Dict{String,Float64}, tp::ThermalParams)
    η = tp.eta_ball
    H_si = get(heat_avg, "H_slide_i", 0.0)
    H_so = get(heat_avg, "H_slide_o", 0.0)
    H_spi = get(heat_avg, "H_spin_i", 0.0)
    H_spo = get(heat_avg, "H_spin_o", 0.0)
    H_drag = get(heat_avg, "H_drag", 0.0)
    H_churn = get(heat_avg, "H_churn", 0.0)

    Q_ir = (1.0 - η) * (H_si + H_spi)
    Q_or = (1.0 - η) * (H_so + H_spo)
    Q_ball = η * (H_si + H_so + H_spi + H_spo)
    Q_oil = H_drag + H_churn

    return (Q_ir, Q_or, Q_ball, Q_oil)
end
