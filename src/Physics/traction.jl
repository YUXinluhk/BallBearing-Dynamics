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
@inline function tehd_traction_coefficient(u_abs, P_mean, b_hz, u_entrain,
    A, B, C, D, Lambda_LSS, beta_temp, effusivity)

    T_ret = typeof(u_abs * P_mean)
    eps_v = T_ret(1e-6)
    u_s_smooth = abs(u_abs) + eps_v

    # 1. Macro isothermal 4-parameter traction (SKF curve)
    mu_iso = (A + B * u_s_smooth) * exp(-C * u_s_smooth) + D

    # Short-circuit if no contact pressure
    (P_mean <= eps_v) && return mu_iso

    # 2. Bair-Winer LSS truncation: branchless, C∞ smooth
    mu_eff = mu_iso / hypot(T_ret(1.0), mu_iso / (Lambda_LSS + eps_v))

    # 3. Archard contact flash temperature
    q_heat = mu_eff * P_mean * u_s_smooth
    dT_flash = (T_ret(1.11) * q_heat * sqrt(abs(b_hz) + eps_v)) /
               (T_ret(effusivity) * sqrt(abs(u_entrain) + eps_v))

    # 4. Nahme thermal softening
    phi_T = T_ret(1.0) / sqrt(T_ret(1.0) + beta_temp * dT_flash)

    # 5. Elastic wave speed limit
    psi_wave = exp(-(u_s_smooth / T_ret(3200.0))^2)

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

const TRACTION_GL7_NODES = (
    0.02544604382862075, 0.12923440720030275, 0.29707742431130140,
    0.5, 0.70292257568869860, 0.87076559279969725, 0.97455395617137925,
)
const TRACTION_GL7_WEIGHTS = (
    0.06474248308443485, 0.13985269574463830, 0.19091502525255945,
    0.20897959183673470, 0.19091502525255945, 0.13985269574463830, 0.06474248308443485,
)
const TRACTION_NTHETA = 12

"""
    integrate_tehd_contact_force(u_roll, u_side, P_mean, a_hz, b_hz, u_entrain,
                                 ω_spin, ω_roll_abs, R_ball, f_geom,
                                 A, B, C, D, Lambda_LSS, beta_temp;
                                 flash_length=sqrt(b_hz))
        → (F_roll, F_side, H_slide, H_spin)

Integrate local TEHD shear stress over the Hertz ellipse and return:
- `F_roll`: signed net traction force along the rolling/sliding direction [N]
- `F_side`: signed net traction force along the lateral creepage direction [N]
- `H_slide`: signed sliding power contribution [W]
- `H_spin`: spin-induced shear power contribution [W]

This is a closer match to Gupta-style local contact mechanics than the older
`κ * Q` lumped traction surrogate.
"""
function integrate_tehd_contact_force(u_roll, u_side, P_mean, a_hz, b_hz, u_entrain,
    ω_spin, ω_roll_abs, R_ball, f_geom,
    A, B, C, D, Lambda_LSS, beta_temp, effusivity;
    flash_length=nothing)

    T_ret = typeof((u_roll + u_side) * P_mean + a_hz * b_hz + ω_spin * R_ball)
    eps_t = T_ret(1e-16)

    (P_mean <= eps_t || a_hz <= eps_t || b_hz <= eps_t || R_ball <= eps_t) &&
        return (zero(T_ret), zero(T_ret), zero(T_ret), zero(T_ret))

    p_max = T_ret(1.5) * P_mean
    conf_factor = (T_ret(2.0) * f_geom - T_ret(1.0)) / (T_ret(2.0) * f_geom)
    dtheta = T_ret(2.0 * π / TRACTION_NTHETA)
    u_roll_eff = sqrt(abs(u_entrain)^2 + T_ret(1e-6))
    flash_len = isnothing(flash_length) ? sqrt(max(b_hz, eps_t)) : sqrt(max(flash_length, eps_t))

    F_roll = zero(T_ret)
    F_side = zero(T_ret)
    H_slide = zero(T_ret)
    H_spin = zero(T_ret)

    for i_u in eachindex(TRACTION_GL7_NODES)
        u_nd = T_ret(TRACTION_GL7_NODES[i_u])
        wu = T_ret(TRACTION_GL7_WEIGHTS[i_u])
        ρ_val = sqrt(T_ret(1.0) - u_nd^2)
        p_local = p_max * u_nd
        τ_lim_eff = sqrt((Lambda_LSS * p_local)^2 + eps_t)

        for i_theta in 1:TRACTION_NTHETA
            θ = (T_ret(i_theta) - T_ret(0.5)) * dtheta
            s_th, c_th = sincos(θ)
            x_c = a_hz * ρ_val * c_th
            y_c = b_hz * ρ_val * s_th

            u_heathcote = ω_roll_abs * (x_c^2 / (T_ret(2.0) * R_ball)) * conf_factor
            v_micro_lat = -ω_spin * y_c
            v_micro_roll = ω_spin * x_c + u_heathcote

            v_lat = u_side + v_micro_lat
            v_roll = u_roll + v_micro_roll
            v_mag = sqrt(v_lat^2 + v_roll^2 + eps_t)

            mu_local_iso = traction_coefficient(v_mag, A, B, C, D)
            τ_iso = mu_local_iso * p_local
            τ_actual = τ_iso / hypot(T_ret(1.0), τ_iso / τ_lim_eff)

            q_local = τ_actual * v_mag
            dT_flash = (T_ret(1.11) * q_local * flash_len) /
                       (T_ret(effusivity) * sqrt(u_roll_eff))
            thermal_reduction = max(T_ret(0.25), T_ret(1.0) / sqrt(T_ret(1.0) + beta_temp * dT_flash))

            τ_final = τ_actual * thermal_reduction
            dF = τ_final * (a_hz * b_hz * u_nd * wu * dtheta)
            dF_lat = dF * (v_lat / v_mag)
            dF_roll = dF * (v_roll / v_mag)

            F_roll += dF_roll
            F_side += dF_lat
            H_slide += dF_roll * u_roll + dF_lat * u_side
            H_spin += dF_lat * v_micro_lat + dF_roll * v_micro_roll
        end
    end

    return (F_roll, F_side, H_slide, H_spin)
end
