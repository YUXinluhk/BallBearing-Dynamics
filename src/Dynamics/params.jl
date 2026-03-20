# =====================================================================
# Dynamics/params.jl — Structured parameter builder for the ODE kernel
#
# Type definitions are in params_types.jl (included before Physics layer).
# This file contains only the build_params constructor function.
# =====================================================================

using Printf
using Statistics: mean

"""
    build_params(geom, mat, lub, trac, cage, config, scales, h_inner, h_outer) → ODEParams

Build structured parameter pack for the ODE kernel (all nondimensional).
"""
function build_params(geom::BearingGeometry, mat::MaterialParams,
    lub::LubricantParams, trac::TractionParams,
    cage::CageGeometry, config::SimulationConfig,
    scales::Scales,
    h_inner::HertzContact, h_outer::HertzContact;
    qs::Union{Nothing,QuasiStaticResult}=nothing)
    s = scales

    # ── Geometry ──
    α_f = alpha_free(geom)
    B_i = (geom.f_i - 0.5) * geom.d
    B_o = (geom.f_o - 0.5) * geom.d
    geom_nd = GeomParams(
        nondim_length(s, geom.d),
        nondim_length(s, geom.d_m),
        geom.f_i, geom.f_o, α_f, geom.n_balls,
        nondim_length(s, D_i(geom)),
        nondim_length(s, D_o(geom)),
        nondim_length(s, (geom.f_i - 0.5) * geom.d),
        nondim_length(s, (geom.f_o - 0.5) * geom.d),
        nondim_length(s, B_i * sin(α_f)),
        nondim_length(s, -B_o * sin(α_f))
    )

    # ── Mass/Inertia ──
    mass_nd = MassParams(
        nondim_mass(s, ball_mass(geom)),
        nondim_inertia(s, ball_inertia(geom)),
        nondim_mass(s, inner_race_mass(geom)),
        nondim_mass(s, cage.cage_mass),
        nondim_inertia(s, cage.cage_inertia_xx),
        nondim_inertia(s, cage.cage_inertia_yy)
    )

    # ── Hertz ──
    hertz_nd = HertzParams(
        composite_modulus(mat, mat) / (s.Q / s.L^2),
        nondim_stiffness(h_inner, s),
        nondim_stiffness(h_outer, s),
        h_inner.a_star, h_inner.b_star,
        h_outer.a_star, h_outer.b_star,
        h_inner.sum_rho * s.L,
        h_outer.sum_rho * s.L,
        h_inner.E2, h_outer.E2,
        config.integrator.eps_contact
    )

    # ── Lubricant ──
    churn = config.churning
    lub_nd = LubParams(
        lub.mu_0, lub.alpha_pv, lub.beta_temp, lub.Lambda_LSS,
        lub.T_0, lub.K_th, lub.rho_lub,
        effective_density(churn.rho_oil, churn.rho_air, churn.fill_fraction),
        churn.mu_oil, mat.effusivity
    )

    # ── Traction ──
    trac_nd = TractionCoeffs(trac.A, trac.B, trac.C, trac.D)

    # ── Per-mode damping ──
    ζ = config.zeta
    Z = geom.n_balls
    m_ball_dim = ball_mass(geom)
    m_ir_dim = inner_race_mass(geom)
    m_cage_dim = cage.cage_mass
    I_ball_dim = ball_inertia(geom)
    r_ball_dim = r_ball(geom)

    if qs !== nothing && any(qs.delta_inner .> 0)
        δ_avg = 0.5 * (mean(qs.delta_inner) + mean(qs.delta_outer))
    else
        δ_avg = 1e-6
    end
    K_hertz = 1.5 * h_inner.Upsilon * sqrt(max(δ_avg, 1e-8))

    c_ball_trans_dim = 2ζ * sqrt(K_hertz * m_ball_dim)
    c_ir_dim = 2ζ * sqrt(Z * K_hertz * m_ir_dim)
    K_cage_proxy = K_hertz * 0.01
    c_cage_dim = 2ζ * sqrt(K_cage_proxy * m_cage_dim)
    c_ball_orbit_dim = 2ζ * sqrt(K_cage_proxy * m_ball_dim)
    K_rot = K_hertz * r_ball_dim^2
    c_ball_spin_dim = 2ζ * sqrt(K_rot * I_ball_dim)
    c_tilt_dim = c_ir_dim * (geom.d_m / 2)^2 * 0.1

    nondim_force_damp(c) = c * s.V / s.Q
    nondim_moment_damp(c) = c * s.W / s.M

    damp_nd = DampingCoeffs(
        nondim_force_damp(c_ball_trans_dim),
        nondim_force_damp(c_ball_orbit_dim),
        nondim_moment_damp(c_ball_spin_dim),
        nondim_force_damp(c_ir_dim),
        nondim_force_damp(c_cage_dim),
        nondim_moment_damp(c_tilt_dim)
    )

    @printf("\n=== Per-Mode Damping (ζ = %.2f) ===\n", ζ)
    @printf("  K_hertz = %.2e N/m  (δ_avg = %.2f µm)\n", K_hertz, δ_avg * 1e6)
    @printf("  c_ball_trans  = %.2f N·s/m  (c★ = %.6f)\n", c_ball_trans_dim, damp_nd.c_ball_trans)
    @printf("  c_ball_orbit  = %.2f N·s/m  (c★ = %.6f)\n", c_ball_orbit_dim, damp_nd.c_ball_orbit)
    @printf("  c_ball_spin   = %.4f N·m·s   (c★ = %.6f)\n", c_ball_spin_dim, damp_nd.c_ball_spin)
    @printf("  c_ir          = %.2f N·s/m  (c★ = %.6f)\n", c_ir_dim, damp_nd.c_ir)
    @printf("  c_cage        = %.2f N·s/m  (c★ = %.6f)\n", c_cage_dim, damp_nd.c_cage)
    @printf("  c_tilt        = %.2f N·s/m  (c★ = %.6f)\n", c_tilt_dim, damp_nd.c_tilt)

    # ── Cage ──
    cage_nd = CageParams(
        nondim_length(s, cage.pocket_clearance),
        cage.stiffness_pocket / (s.Q / s.L),
        cage.mu_pocket,
        cage.stiffness_pilot / (s.Q / s.L),
        cage.mu_pilot,
        nondim_length(s, cage.pilot_clearance),
        cage.pilot_is_inner,
        nondim_length(s, cage.cage_inner_radius),
        nondim_length(s, cage.cage_outer_radius),
        nondim_length(s, churn.cage_web),
        cage.c_damping / (s.Q * s.T / s.L)
    )

    # ── Loads ──
    ω_ir_dim = config.inner_race_speed
    γ = geom.d * cos(α_f) / geom.d_m
    ω_cage_dim = 0.5 * ω_ir_dim * (1 - γ)

    load_nd = LoadParams(
        nondim_angvel(s, config.inner_race_speed),
        nondim_angvel(s, config.outer_race_speed),
        nondim_angvel(s, ω_cage_dim),
        nondim_force(s, config.F_axial),
        nondim_force(s, config.F_radial),
        nondim_time(s, config.t_ramp_end),
        config.mu_spin, config.zeta
    )

    # ── Scales ──
    scale_nd = ScaleFactors(s.V, s.W, s.L, s.Q, s.T)

    # ── Thermal ──
    th = config.thermal
    th_accel = th.th_accel
    C_ir = (isnan(th.C_ir) ? inner_race_mass(geom) * th.c_steel : th.C_ir) / th_accel
    C_or = (isnan(th.C_or) ? inner_race_mass(geom) * 1.5 * th.c_steel : th.C_or) / th_accel
    C_ball = (isnan(th.C_ball) ? ball_mass(geom) * Z * th.c_steel : th.C_ball) / th_accel
    C_oil = (isnan(th.C_oil) ? ball_mass(geom) * Z * 10.0 * (lub.rho_lub / geom.rho_ball) * lub.c_p : th.C_oil) / th_accel

    nondim_mcp(C) = C / (s.Q * s.L)
    nondim_G(G) = G / (s.Q * s.L / s.T)

    T_oil_in = isnan(th.T_oil_inlet) ? th.T_init : th.T_oil_inlet
    if th.oil_flow_rate > 0.0
        V_dot_m3s = th.oil_flow_rate / 60.0 / 1e6
        m_dot = lub.rho_lub * V_dot_m3s
        mdot_cp = m_dot * lub.c_p
    else
        mdot_cp = 0.0
    end

    thermal_nd = ThermalNDParams(
        nondim_mcp(C_ir), nondim_mcp(C_or), nondim_mcp(C_ball), nondim_mcp(C_oil),
        nondim_G(th.G_ir_ball), nondim_G(th.G_or_ball),
        nondim_G(th.G_ball_oil), nondim_G(th.G_or_amb), nondim_G(th.G_oil_amb),
        th.T_ambient, th.CTE, th.T_init,
        geom.rho_ball, mat.E,
        nondim_G(mdot_cp), T_oil_in
    )

    return ODEParams(geom_nd, mass_nd, hertz_nd, lub_nd, trac_nd, damp_nd,
                     cage_nd, load_nd, scale_nd, thermal_nd)
end

# ── Legacy compat: keep N_PARAMS symbol for postprocess backward compat ──
const N_PARAMS = 93
