# =====================================================================
# Dynamics/params.jl — Flat parameter vector for the ODE kernel
#
# Port of dynamics.py::_build_numba_params.
# Named constants for indexing into the parameter array.
# =====================================================================

using Printf
using Statistics: mean

# ── Parameter indices (1-indexed) ─────────────────────────────────────

const P_D = 1   # ball diameter
const P_DM = 2   # pitch diameter
const P_FI = 3   # inner conformity
const P_FO = 4   # outer conformity
const P_ALPHA0 = 5   # free contact angle
const P_NBALL = 6   # number of balls
const P_MBALL = 7   # ball mass [nondim]
const P_JBALL = 8   # ball inertia [nondim]
const P_MIR = 9   # inner race mass [nondim]
const P_E_PRIME = 10  # composite modulus [nondim]
const P_YI = 11  # inner Hertz Υ★
const P_YO = 12  # outer Hertz Υ★
const P_AI = 13  # inner a★
const P_BO = 14  # outer b★
const P_BI = 15  # inner b★
const P_AO = 16  # outer a★
const P_SRI = 17  # Σρ inner [nondim]
const P_SRO = 18  # Σρ outer [nondim]
const P_DI = 19  # D_inner [nondim]
const P_DO = 20  # D_outer [nondim]
const P_MU0 = 21  # μ₀ viscosity [nondim]
const P_ALPHA_PV = 22  # α_pv
const P_BETA_TEMP = 23  # β_temp
const P_LAMBDA_LSS = 24  # Λ_LSS
const P_T0 = 25  # T₀
const P_KTH = 26  # K_th
const P_RHO_LUB = 27  # ρ_lub
const P_TRAC_A = 28  # traction A
const P_TRAC_B = 29  # traction B
const P_TRAC_C = 30  # traction C
const P_TRAC_D = 31  # traction D
const P_MU_SPIN = 32  # μ_spin
const P_C_STRUCT = 33  # structural damping (legacy, may be overridden by per-mode)
const P_ZETA = 34  # damping ratio
const P_OMEGA_IR = 35  # inner race speed [nondim]
const P_FA = 36  # axial load [nondim]
const P_FR = 37  # radial load [nondim]
const P_T_RAMP = 38  # ramp time [nondim]
const P_RHO_EFF = 39  # effective density [nondim]
const P_MU_OIL = 40  # oil viscosity [nondim]
const P_CAGE_WEB = 41  # cage web [nondim]
const P_POCKET_R = 42  # pocket radius [nondim]
const P_POCKET_CLR = 43  # pocket clearance [nondim]
const P_CAGE_MASS = 44  # cage mass [nondim]
const P_CAGE_IXX = 45  # cage Ixx [nondim]
const P_CAGE_IYY = 46  # cage Iyy [nondim]
const P_K_POCKET = 47  # pocket stiffness [nondim]
const P_K_PILOT = 48  # pilot stiffness [nondim]
const P_MU_POCKET = 49  # pocket friction
const P_MU_PILOT = 50  # pilot friction
const P_C_CAGE = 51  # cage damping
const P_PILOT_CLR = 52  # pilot clearance [nondim]
const P_PILOT_IR = 53  # pilot is inner (1/0)
const P_CAGE_IR = 54  # cage inner radius [nondim]
const P_CAGE_OR = 55  # cage outer radius [nondim]
const P_DRB_I = 56  # distance race-ball center inner
const P_DRB_O = 57  # distance race-ball center outer
const P_DELTA_R_TH = 58  # thermal radial expansion [nondim]
const P_X_GI0 = 59  # inner groove center axial offset ★
const P_X_GO0 = 60  # outer groove center axial offset ★
# Per-mode damping (matching Python dynamics.py per-mode model)
const P_C_BALL_TRANS = 61  # ball axial/radial damping ★
const P_C_BALL_ORBIT = 62  # ball orbital θ damping ★
const P_C_BALL_SPIN = 63  # ball spin ω damping ★
const P_C_IR_DAMP = 64  # inner race translation damping ★
const P_C_CAGE_DAMP = 65  # cage translation/rotation damping ★
const P_C_TILT = 66  # inner race tilt damping ★
const P_OMEGA_CAGE = 67  # kinematic cage speed ★
const P_V_SCALE = 68     # velocity scale V [m/s] (for traction/drag dimensionalization)
const P_W_SCALE = 69     # angular velocity scale W [rad/s] (for spin moment dimensionalization)
const P_L_SCALE = 70     # length scale L [m] (for contact geometry dimensionalization)
const P_Q_SCALE = 71     # force scale Q [N] (for force dimensionalization)
const P_E2_I = 72        # E(m) inner contact (precomputed elliptic integral)
const P_E2_O = 73        # E(m) outer contact (precomputed elliptic integral)
const N_PARAMS = 73

"""
    build_params(geom, mat, lub, trac, cage, config, scales, h_inner, h_outer) → Vector{Float64}

Build flat parameter vector for the ODE kernel (all non-dimensional).
"""
function build_params(geom::BearingGeometry, mat::MaterialParams,
    lub::LubricantParams, trac::TractionParams,
    cage::CageGeometry, config::SimulationConfig,
    scales::Scales,
    h_inner::HertzContact, h_outer::HertzContact;
    qs::Union{Nothing,QuasiStaticResult}=nothing)
    s = scales
    p = zeros(N_PARAMS)

    p[P_D] = nondim_length(s, geom.d)
    p[P_DM] = nondim_length(s, geom.d_m)
    p[P_FI] = geom.f_i
    p[P_FO] = geom.f_o
    p[P_ALPHA0] = geom.alpha_0
    p[P_NBALL] = Float64(geom.n_balls)
    p[P_MBALL] = nondim_mass(s, ball_mass(geom))
    p[P_JBALL] = nondim_inertia(s, ball_inertia(geom))
    p[P_MIR] = nondim_mass(s, inner_race_mass(geom))
    p[P_E_PRIME] = composite_modulus(mat, mat) / (s.Q / s.L^2)
    p[P_YI] = nondim_stiffness(h_inner, s)
    p[P_YO] = nondim_stiffness(h_outer, s)
    p[P_AI] = h_inner.a_star
    p[P_BI] = h_inner.b_star
    p[P_AO] = h_outer.a_star
    p[P_BO] = h_outer.b_star
    p[P_SRI] = h_inner.sum_rho * s.L
    p[P_SRO] = h_outer.sum_rho * s.L
    p[P_DI] = nondim_length(s, D_i(geom))
    p[P_DO] = nondim_length(s, D_o(geom))
    p[P_MU0] = lub.mu_0
    p[P_ALPHA_PV] = lub.alpha_pv
    p[P_BETA_TEMP] = lub.beta_temp
    p[P_LAMBDA_LSS] = lub.Lambda_LSS
    p[P_T0] = lub.T_0
    p[P_KTH] = lub.K_th
    p[P_RHO_LUB] = lub.rho_lub
    p[P_TRAC_A] = trac.A
    p[P_TRAC_B] = trac.B
    p[P_TRAC_C] = trac.C
    p[P_TRAC_D] = trac.D
    p[P_MU_SPIN] = config.mu_spin
    p[P_C_STRUCT] = config.c_structural / (s.Q * s.T / s.L)  # damping nondim
    p[P_ZETA] = config.zeta
    p[P_OMEGA_IR] = nondim_angvel(s, config.inner_race_speed)
    p[P_FA] = nondim_force(s, config.F_axial)
    p[P_FR] = nondim_force(s, config.F_radial)
    p[P_T_RAMP] = nondim_time(s, config.t_ramp_end)
    churn = config.churning
    p[P_RHO_EFF] = effective_density(churn.rho_oil, churn.rho_air, churn.fill_fraction)
    p[P_MU_OIL] = churn.mu_oil
    p[P_CAGE_WEB] = nondim_length(s, churn.cage_web)
    p[P_POCKET_R] = nondim_length(s, cage.pocket_radius)
    p[P_POCKET_CLR] = nondim_length(s, cage.pocket_clearance)
    p[P_CAGE_MASS] = nondim_mass(s, cage.cage_mass)
    p[P_CAGE_IXX] = nondim_inertia(s, cage.cage_inertia_xx)
    p[P_CAGE_IYY] = nondim_inertia(s, cage.cage_inertia_yy)
    p[P_K_POCKET] = cage.stiffness_pocket / (s.Q / s.L)
    p[P_K_PILOT] = cage.stiffness_pilot / (s.Q / s.L)
    p[P_MU_POCKET] = cage.mu_pocket
    p[P_MU_PILOT] = cage.mu_pilot
    p[P_C_CAGE] = cage.c_damping / (s.Q * s.T / s.L)
    p[P_PILOT_CLR] = nondim_length(s, cage.pilot_clearance)
    p[P_PILOT_IR] = cage.pilot_is_inner ? 1.0 : 0.0
    p[P_CAGE_IR] = nondim_length(s, cage.cage_inner_radius)
    p[P_CAGE_OR] = nondim_length(s, cage.cage_outer_radius)
    p[P_DRB_I] = nondim_length(s, (geom.f_i - 0.5) * geom.d)
    p[P_DRB_O] = nondim_length(s, (geom.f_o - 0.5) * geom.d)
    p[P_DELTA_R_TH] = nondim_length(s, config.delta_r_thermal)
    # Groove center axial offsets (Harris §3.5): B·sin(α₀)
    B_i = (geom.f_i - 0.5) * geom.d
    B_o = (geom.f_o - 0.5) * geom.d
    p[P_X_GI0] = nondim_length(s, B_i * sin(geom.alpha_0))  # positive
    p[P_X_GO0] = nondim_length(s, -B_o * sin(geom.alpha_0))  # negative

    # ── Per-mode damping (matching Python dynamics.py L900-923) ──
    ζ = config.zeta
    Z = geom.n_balls
    m_ball_dim = ball_mass(geom)
    m_ir_dim = inner_race_mass(geom)
    m_cage_dim = cage.cage_mass
    I_ball_dim = ball_inertia(geom)
    r_ball_dim = r_ball(geom)

    # Average Hertz stiffness from QS penetration
    if qs !== nothing && any(qs.delta_inner .> 0)
        δ_avg = 0.5 * (mean(qs.delta_inner) + mean(qs.delta_outer))
    else
        δ_avg = 1e-6  # fallback
    end
    K_hertz = 1.5 * h_inner.Upsilon * sqrt(max(δ_avg, 1e-8))

    c_ball_trans_dim = 2ζ * sqrt(K_hertz * m_ball_dim)
    c_ir_dim = 2ζ * sqrt(Z * K_hertz * m_ir_dim)
    K_cage_proxy = K_hertz * 0.01
    c_cage_dim = 2ζ * sqrt(K_cage_proxy * m_cage_dim)
    c_ball_orbit_dim = 2ζ * sqrt(K_cage_proxy * m_ball_dim)
    K_rot = K_hertz * r_ball_dim^2
    c_ball_spin_dim = 2ζ * sqrt(K_rot * I_ball_dim)
    c_tilt_dim = max(c_ir_dim, 2ζ * sqrt(K_hertz * m_ir_dim * (geom.d_m / 2)^2))

    # Nondimensionalize: c★ = c · V / Q₀  (force-damping), or c · W / M₀ (moment-damping)
    nondim_damp(c) = c * s.V / s.Q
    p[P_C_BALL_TRANS] = nondim_damp(c_ball_trans_dim)
    p[P_C_BALL_ORBIT] = nondim_damp(c_ball_orbit_dim)
    p[P_C_BALL_SPIN] = nondim_damp(c_ball_spin_dim)
    p[P_C_IR_DAMP] = nondim_damp(c_ir_dim)
    p[P_C_CAGE_DAMP] = nondim_damp(c_cage_dim)
    p[P_C_TILT] = max(nondim_damp(c_tilt_dim), 1.0)

    # Kinematic cage speed for relative orbital damping
    ω_ir_dim = config.inner_race_speed
    γ = geom.d * cos(geom.alpha_0) / geom.d_m
    ω_cage_dim = 0.5 * ω_ir_dim * (1 - γ)
    p[P_OMEGA_CAGE] = nondim_angvel(s, ω_cage_dim)

    @printf("\n=== Per-Mode Damping (ζ = %.2f) ===\n", ζ)
    @printf("  K_hertz = %.2e N/m  (δ_avg = %.2f µm)\n", K_hertz, δ_avg * 1e6)
    @printf("  c_ball_trans  = %.2f N·s/m  (c★ = %.6f)\n", c_ball_trans_dim, p[P_C_BALL_TRANS])
    @printf("  c_ball_orbit  = %.2f N·s/m  (c★ = %.6f)\n", c_ball_orbit_dim, p[P_C_BALL_ORBIT])
    @printf("  c_ball_spin   = %.4f N·m·s   (c★ = %.6f)\n", c_ball_spin_dim, p[P_C_BALL_SPIN])
    @printf("  c_ir          = %.2f N·s/m  (c★ = %.6f)\n", c_ir_dim, p[P_C_IR_DAMP])
    @printf("  c_cage        = %.2f N·s/m  (c★ = %.6f)\n", c_cage_dim, p[P_C_CAGE_DAMP])
    @printf("  c_tilt        = %.2f N·s/m  (c★ = %.6f)\n", c_tilt_dim, p[P_C_TILT])

    # Dimensional scales for kernel (traction/drag/spin use dimensional inputs)
    p[P_V_SCALE] = s.V
    p[P_W_SCALE] = s.W
    p[P_L_SCALE] = s.L
    p[P_Q_SCALE] = s.Q

    # Precomputed elliptic integrals E(m) for spin moment (avoids hot-loop ellipE calls)
    p[P_E2_I] = h_inner.E2
    p[P_E2_O] = h_outer.E2

    return p
end
