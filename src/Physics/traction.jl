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

const TRACTION_GL5_NODES = (
    0.0469100770306680, 0.23076534494715845, 0.5,
    0.76923465505284155, 0.9530899229693320,
)
const TRACTION_GL5_WEIGHTS = (
    0.11846344252809455, 0.23931433524968325, 0.28444444444444444,
    0.23931433524968325, 0.11846344252809455,
)
const TRACTION_NTHETA = 8

const TRACTION_GL5_RHO = (
    0.21658734272895963, 0.48038041690628863, 0.7071067811865476,
    0.8770602345636413, 0.9762632447036668,
)
const TRACTION_SIN_THETA = (
    0.7071067811865475, 1.0, 0.7071067811865476, 1.2246467991473532e-16,
    -0.7071067811865475, -1.0, -0.7071067811865476, -2.4492935982947064e-16,
)
const TRACTION_COS_THETA = (
    0.7071067811865476, 6.123233995736766e-17, -0.7071067811865475, -1.0,
    -0.7071067811865477, -1.8369701987210297e-16, 0.7071067811865474, 1.0,
)

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
- `M_spin`: net spin friction torque [N·m]

This is a closer match to Gupta-style local contact mechanics than the older
`κ * Q` lumped traction surrogate.
"""
function integrate_tehd_contact_force(u_roll, u_side, P_mean, a_hz, b_hz, u_entrain,
    ω_spin, ω_roll_abs, R_ball, f_geom,
    A, B, C, D, Lambda_LSS, beta_temp, effusivity;
    flash_length=nothing)

    T_ret = typeof((u_roll + u_side) * P_mean + a_hz * b_hz + ω_spin * R_ball)
    eps_num = T_ret(1e-12)

    (P_mean <= eps_num || a_hz <= eps_num || b_hz <= eps_num || R_ball <= eps_num) &&
        return (zero(T_ret), zero(T_ret), zero(T_ret), zero(T_ret), zero(T_ret))

    p_max = T_ret(1.5) * P_mean
    conf_factor = (T_ret(2.0) * f_geom - T_ret(1.0)) / (T_ret(2.0) * f_geom)
    dtheta = T_ret(2.0 * π / TRACTION_NTHETA)
    
    u_roll_eff = sqrt(u_entrain^2 + T_ret(1e-6))
    flash_len = isnothing(flash_length) ? sqrt(b_hz^2 + eps_num) : sqrt(flash_length^2 + eps_num)

    F_roll = zero(T_ret)
    F_side = zero(T_ret)
    H_slide = zero(T_ret)
    H_spin = zero(T_ret)
    M_spin = zero(T_ret)

    heathcote_coeff = ω_roll_abs * conf_factor / (T_ret(2.0) * R_ball)
    flash_coeff = T_ret(1.11) * flash_len / (T_ret(effusivity) * sqrt(u_roll_eff))

    # 循环顺序重构：允许 CPU 利用 @fastmath 展开内循环
    @fastmath @inbounds for i_u in 1:5
        u_nd = T_ret(TRACTION_GL5_NODES[i_u])
        wu = T_ret(TRACTION_GL5_WEIGHTS[i_u])
        ρ_val = T_ret(TRACTION_GL5_RHO[i_u])
        p_local = p_max * sqrt(max(zero(T_ret), T_ret(1.0) - u_nd))
        τ_lim_eff_sq = (Lambda_LSS * p_local)^2 + eps_num

        @inbounds for i_theta in 1:8
            s_th = T_ret(TRACTION_SIN_THETA[i_theta])
            c_th = T_ret(TRACTION_COS_THETA[i_theta])
            
            x_c = a_hz * c_th * ρ_val
            y_c = b_hz * s_th * ρ_val

            u_heathcote = heathcote_coeff * x_c^2
            v_micro_lat = -ω_spin * y_c
            v_micro_roll = ω_spin * x_c + u_heathcote

            v_lat = u_side + v_micro_lat
            v_roll = u_roll + v_micro_roll
            
            # v_mag 自带 eps_num，完美锁死隐式求导
            v_mag = sqrt(v_lat^2 + v_roll^2 + eps_num)

            mu_local_iso = (A + B * v_mag) * exp(-C * v_mag) + D + T_ret(1e-6) * v_mag
            τ_iso = mu_local_iso * p_local
            
            # 解析打通 Bair-Winer 平滑
            τ_actual = τ_iso / sqrt(T_ret(1.0) + τ_iso^2 / τ_lim_eff_sq)

            q_local = τ_actual * v_mag
            dT_flash = flash_coeff * q_local
            
            thermal_reduction = T_ret(1.0) / sqrt(T_ret(1.0) + beta_temp * dT_flash)

            τ_final = τ_actual * thermal_reduction
            # Jacobian: polar u=ρ² ⟹ ρ dρ = ½ du, so dA = a·b·½·du·dθ
            dF = τ_final * (a_hz * b_hz * 0.5 * wu * dtheta)
            
            dF_lat = dF * (v_lat / v_mag)
            dF_roll = dF * (v_roll / v_mag)

            F_roll += dF_roll
            F_side += dF_lat
            H_slide += dF_roll * u_roll + dF_lat * u_side
            H_spin += dF_lat * v_micro_lat + dF_roll * v_micro_roll
            M_spin += x_c * dF_roll - y_c * dF_lat
        end
    end

    return (F_roll, F_side, H_slide, H_spin, M_spin)
end
