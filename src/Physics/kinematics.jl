# =====================================================================
# Physics/kinematics.jl — Core TEHD Kinematics & EHL Forces
#
# Shared purely mathematical logic for both the ODE Kernel and Field Output.
# Decouples geometric tracking, velocities, Hertz, and EHL from the solver.
# =====================================================================

using LinearAlgebra
using StaticArrays

"""
    KinematicContext

Stores all system-level (non-ball-specific) scalars and vectors needed
for calculating ball contact kinematics. Precomputed once per time step.
"""
struct KinematicContext{T}
    u::AbstractVector{T}
    p::ODEParams
    ramp::Float64
    ω_ir::Float64
    ω_or::Float64
    
    # Scales
    V_scale::Float64
    W_scale::Float64
    L_scale::Float64
    Q_scale::Float64
    E_prime_dim::Float64
    ω_ir_dim::Float64
    d_m_star::Float64
    
    # Ring states
    x_ir::T; y_ir::T; z_ir::T; γy_ir::T; γz_ir::T
    ẋ_ir::T; ẏ_ir::T; ż_ir::T; γ̇y_ir::T; γ̇z_ir::T
    x_cg::T; y_cg::T; z_cg::T; θ_cg::T
    ẋ_cg::T; ẏ_cg::T; ż_cg::T; θ̇_cg::T
    
    # Thermal states
    T_ir_node::T; T_or_node::T; T_ball_node::T; T_oil_node::T
    
    # Derived EHL & Thermo properties
    μ_oil_dyn::T
    ν_cst_064::T
    G_val::Float64
    eps_c::Float64
    entrain_factor::T
end

@inline function build_kinematic_context(u::AbstractVector, p::ODEParams, t::Float64)
    T_u = eltype(u)
    Z = p.geom.Z
    
    t_ramp = p.load.t_ramp
    ramp = t_ramp > 0 ? 0.5 * (1.0 - cos(π * min(t / t_ramp, 1.0))) : 1.0
    ω_ir = p.load.omega_ir * ramp
    ω_or = p.load.omega_or * ramp
    
    V_scale = p.scale.V
    W_scale = p.scale.W
    L_scale = p.scale.L
    Q_scale = p.scale.Q
    d_m_star = p.geom.d_m
    E_prime_dim = p.hertz.E_prime * Q_scale / L_scale^2
    
    ir_p = ir_pos_view(u, Z)
    ir_v = ir_vel_view(u, Z)
    x_ir, y_ir, z_ir, γy_ir, γz_ir = Tuple(ir_p)
    ẋ_ir, ẏ_ir, ż_ir, γ̇y_ir, γ̇z_ir = Tuple(ir_v)
    
    cg_p = cage_pos_view(u, Z)
    cg_v = cage_vel_view(u, Z)
    x_cg, y_cg, z_cg, θ_cg = Tuple(cg_p)
    ẋ_cg, ẏ_cg, ż_cg, θ̇_cg = Tuple(cg_v)
    
    tv = thermal_view(u, Z)
    T_ir_node, T_or_node, T_ball_node, T_oil_node = Tuple(tv)
    
    T_ref_oil = p.lub.T_0
    μ_oil_dyn = p.lub.mu_oil * exp(clamp(-p.lub.beta_temp * (T_oil_node - T_ref_oil), -30.0, 30.0))
    ν_cst = p.lub.rho_eff > 0.0 ? (μ_oil_dyn / p.lub.rho_eff) * 1e6 : T_u(10.0)
    ν_cst_064 = ν_cst^0.64
    
    G_val = p.lub.alpha_pv * p.hertz.E_prime
    ω_ir_dim = ω_ir * W_scale
    
    v_roll_nominal_sq = (ω_ir_dim * d_m_star * L_scale * 0.5)^2
    entrain_factor = v_roll_nominal_sq / (v_roll_nominal_sq + 0.01)
    
    return KinematicContext{T_u}(
        u, p, ramp, ω_ir, ω_or,
        V_scale, W_scale, L_scale, Q_scale, E_prime_dim, ω_ir_dim, d_m_star,
        x_ir, y_ir, z_ir, γy_ir, γz_ir,
        ẋ_ir, ẏ_ir, ż_ir, γ̇y_ir, γ̇z_ir,
        x_cg, y_cg, z_cg, θ_cg,
        ẋ_cg, ẏ_cg, ż_cg, θ̇_cg,
        T_ir_node, T_or_node, T_ball_node, T_oil_node,
        μ_oil_dyn, ν_cst_064, G_val, p.hertz.eps_contact, entrain_factor
    )
end

@inline function compute_ball_tehd_kinematics(ctx::KinematicContext{T_u}, j::Int) where {T_u}
    u = ctx.u
    p = ctx.p
    Z = p.geom.Z
    
    bp = ball_pos_view(u, j, Z)
    bv = ball_vel_view(u, j, Z)
    x_b, r_b, θ_b = bp[1], bp[2], bp[3]
    ẋ_b, ṙ_b, θ̇_b = bv[1], bv[2], bv[3]
    
    ω_body = ball_omega(u, j, Z)
    sθ, cθ = sincos(θ_b)
    
    D_star = p.geom.D
    D_i_star = p.geom.D_i
    D_o_star = p.geom.D_o
    x_gi0_star = p.geom.x_gi0
    x_go0_star = p.geom.x_go0
    drb_i = p.geom.drb_i
    drb_o = p.geom.drb_o
    
    # Phase 4b: Thermal Expansion
    R_ball_k_nom = D_star / 2.0
    δr_ball = ctx.ramp * p.thermal.CTE * (ctx.T_ball_node - p.thermal.T_ref) * R_ball_k_nom
    R_ball_k = R_ball_k_nom + δr_ball
    
    δr_ir = ctx.ramp * p.thermal.CTE * (ctx.T_ir_node - p.thermal.T_ref) * (D_i_star / 2.0)
    δr_or = ctx.ramp * p.thermal.CTE * (ctx.T_or_node - p.thermal.T_ref) * (D_o_star / 2.0)
    R_ir_dim = D_i_star / 2.0 * ctx.L_scale
    ω_ir_dim_local = ctx.ω_ir * ctx.W_scale
    δr_centrifugal = p.thermal.rho_race * ω_ir_dim_local^2 * R_ir_dim^3 / (p.thermal.E_young + 1e-30) / ctx.L_scale
    
    R_pitch_i = D_i_star / 2.0 + δr_ir + δr_centrifugal
    tilt_proj = cθ * ctx.γy_ir + sθ * ctx.γz_ir
    
    x_gc_i = x_gi0_star + ctx.x_ir + R_pitch_i * tilt_proj
    r_gc_i = R_pitch_i + ctx.y_ir * (-sθ) + ctx.z_ir * cθ - x_gi0_star * tilt_proj
    x_gc_o = x_go0_star
    r_gc_o = D_o_star / 2.0 + δr_or
    
    dx_i = x_gc_i - x_b
    dr_i = r_gc_i - r_b
    dx_o = x_gc_o - x_b
    dr_o = r_gc_o - r_b
    
    a2_i_state = zero(T_u)
    a2_o_state = zero(T_u)
    
    L_i_plane = hypot(dx_i, dr_i)
    L_o_plane = hypot(dx_o, dr_o)
    dθ_contact_i = L_i_plane * a2_i_state
    dθ_contact_o = L_o_plane * a2_o_state
    
    n_i_loc, t_lat_i_loc, t_roll_i_loc = contact_basis_from_displacements(dx_i, dr_i, dθ_contact_i)
    n_o_loc, t_lat_o_loc, t_roll_o_loc = contact_basis_from_displacements(dx_o, dr_o, dθ_contact_o)
    
    L_i = sqrt(dx_i^2 + dr_i^2 + dθ_contact_i^2 + 1e-30)
    L_o = sqrt(dx_o^2 + dr_o^2 + dθ_contact_o^2 + 1e-30)
    
    nfi_x, nfi_r, nfi_θ = n_i_loc[1], n_i_loc[2], n_i_loc[3]
    nfo_x, nfo_r, nfo_θ = n_o_loc[1], n_o_loc[2], n_o_loc[3]
    
    δ_i_sm = smooth_hertz_delta(L_i - drb_i; ε=ctx.eps_c)
    δ_o_sm = smooth_hertz_delta(L_o - drb_o; ε=ctx.eps_c)
    Q_i = smooth_hertz_load(L_i - drb_i, p.hertz.Y_i; ε=ctx.eps_c)
    Q_o = smooth_hertz_load(L_o - drb_o, p.hertz.Y_o; ε=ctx.eps_c)
    Q_i_dim = Q_i * ctx.Q_scale
    Q_o_dim = Q_o * ctx.Q_scale
    
    c_i = cbrt(3.0 * Q_i / (2.0 * p.hertz.E_prime * p.hertz.sr_i + 1e-30) + 1e-24)
    a_i_dim, b_i_dim = p.hertz.a_i * c_i * ctx.L_scale, p.hertz.b_i * c_i * ctx.L_scale
    c_o = cbrt(3.0 * Q_o / (2.0 * p.hertz.E_prime * p.hertz.sr_o + 1e-30) + 1e-24)
    a_o_dim, b_o_dim = p.hertz.a_o * c_o * ctx.L_scale, p.hertz.b_o * c_o * ctx.L_scale
    
    r_contact_i = r_b - R_ball_k * nfi_r
    x_contact_i = x_b - R_ball_k * nfi_x
    r_contact_o = r_b - R_ball_k * nfo_r
    x_contact_o = x_b - R_ball_k * nfo_x
    
    v_ir_θ = ctx.ω_ir * (r_contact_i + ctx.y_ir * sθ - ctx.z_ir * cθ) + ctx.γ̇y_ir * (x_contact_i - ctx.x_ir) * sθ - ctx.γ̇z_ir * (x_contact_i - ctx.x_ir) * cθ - ctx.ẏ_ir * cθ - ctx.ż_ir * sθ
    v_ball_θ = r_b * θ̇_b
    v_or_θ = ctx.ω_or * r_contact_o
    
    # ── 1. Exact 3D Contact-Basis Kinematics (Cross-Track) ──
    X_rel_i = x_contact_i - ctx.x_ir
    Y_rel_i = -r_contact_i * sθ - ctx.y_ir
    Z_rel_i = r_contact_i * cθ - ctx.z_ir
    
    v_ir_X_g = ctx.ẋ_ir + ctx.γ̇y_ir * Z_rel_i - ctx.γ̇z_ir * Y_rel_i
    v_ir_Y_g = ctx.ẏ_ir + ctx.γ̇z_ir * X_rel_i - ctx.ω_ir * Z_rel_i
    v_ir_Z_g = ctx.ż_ir + ctx.ω_ir * Y_rel_i - ctx.γ̇y_ir * X_rel_i
    
    v_ir_x = v_ir_X_g
    v_ir_r = -v_ir_Y_g * sθ + v_ir_Z_g * cθ
    v_ir_θ_cyl = -v_ir_Y_g * cθ - v_ir_Z_g * sθ
    
    ω_ir_r_cyl = -ctx.γ̇y_ir * sθ + ctx.γ̇z_ir * cθ
    ω_ir_θ_cyl = -ctx.γ̇y_ir * cθ - ctx.γ̇z_ir * sθ
    
    v_cp_i_dim = r_contact_i * θ̇_b * ctx.V_scale
    v_cp_o_dim = r_contact_o * θ̇_b * ctx.V_scale
    
    v_ir_θ_dim = v_ir_θ_cyl * ctx.V_scale
    v_or_θ_dim = v_or_θ * ctx.V_scale
    R_ball_dim = R_ball_k * ctx.L_scale
    
    v_ball_center_loc = SVector{3,T_u}(ẋ_b * ctx.V_scale, ṙ_b * ctx.V_scale, v_ball_θ * ctx.V_scale)
    v_ir_loc = SVector{3,T_u}(v_ir_x * ctx.V_scale, v_ir_r * ctx.V_scale, v_ir_θ_dim)
    v_or_loc = SVector{3,T_u}(zero(T_u), zero(T_u), v_or_θ_dim)
    
    ω_ball_loc = SVector{3,T_u}(ω_body[1] * ctx.W_scale, ω_body[2] * ctx.W_scale, ω_body[3] * ctx.W_scale)
    ω_ir_loc = SVector{3,T_u}(ctx.ω_ir_dim, ω_ir_r_cyl * ctx.W_scale, ω_ir_θ_cyl * ctx.W_scale)
    ω_or_loc = SVector{3,T_u}(ctx.ω_or * ctx.W_scale, zero(T_u), zero(T_u))
    
    ω_ball_lat_i = dot(ω_ball_loc, t_lat_i_loc)
    ω_ball_lat_o = dot(ω_ball_loc, t_lat_o_loc)
    ω_ball_roll_i = dot(ω_ball_loc, t_roll_i_loc)
    ω_ball_roll_o = dot(ω_ball_loc, t_roll_o_loc)
    
    v_spin_i_loc = R_ball_dim * (ω_ball_lat_i * t_roll_i_loc - ω_ball_roll_i * t_lat_i_loc)
    v_spin_o_loc = R_ball_dim * (ω_ball_lat_o * t_roll_o_loc - ω_ball_roll_o * t_lat_o_loc)
    v_ball_contact_i_loc = v_ball_center_loc + v_spin_i_loc
    v_ball_contact_o_loc = v_ball_center_loc + v_spin_o_loc
    
    rel_v_i_loc = v_ir_loc - v_ball_contact_i_loc
    rel_v_o_loc = v_or_loc - v_ball_contact_o_loc
    u_side_i_raw = dot(rel_v_i_loc, t_lat_i_loc)
    u_side_o_raw = dot(rel_v_o_loc, t_lat_o_loc)
    u_roll_i_raw = dot(rel_v_i_loc, t_roll_i_loc)
    u_roll_o_raw = dot(rel_v_o_loc, t_roll_o_loc)
    
    v_cp_i_transport = v_cp_i_dim * t_roll_i_loc[3]
    v_cp_o_transport = v_cp_o_dim * t_roll_o_loc[3]
    v_ir_roll = dot(v_ir_loc, t_roll_i_loc)
    v_or_roll = dot(v_or_loc, t_roll_o_loc)
    v_ball_roll_i = dot(v_ball_contact_i_loc, t_roll_i_loc)
    v_ball_roll_o = dot(v_ball_contact_o_loc, t_roll_o_loc)
    
    u_entrain_i = hypot(0.5 * (v_ir_roll + v_ball_roll_i) - v_cp_i_transport, T_u(1e-12))
    u_entrain_o = hypot(0.5 * (v_or_roll + v_ball_roll_o) - v_cp_o_transport, T_u(1e-12))
    
    u_roll_i_dim = u_roll_i_raw
    u_roll_o_dim = u_roll_o_raw
    u_side_i_dim = u_side_i_raw
    u_side_o_dim = u_side_o_raw
    
    ω_spin_i = (ω_ir_loc[1] - ω_ball_loc[1]) * nfi_x + (ω_ir_loc[2] - ω_ball_loc[2]) * nfi_r + (ω_ir_loc[3] - ω_ball_loc[3]) * nfi_θ
    ω_spin_o = (ω_or_loc[1] - ω_ball_loc[1]) * nfo_x + (ω_or_loc[2] - ω_ball_loc[2]) * nfo_r + (ω_or_loc[3] - ω_ball_loc[3]) * nfo_θ
    ω_roll_abs_i = sqrt(ω_ball_lat_i^2 + (v_ball_θ / (R_ball_k + 1e-30) * ctx.W_scale)^2 + 1e-16)
    ω_roll_abs_o = sqrt(ω_ball_lat_o^2 + (v_ball_θ / (R_ball_k + 1e-30) * ctx.W_scale)^2 + 1e-16)
    
    # ── 2. 2D TEHD Contact Integration ──
    P_mean_i = Q_i_dim / (π * a_i_dim * b_i_dim + 1e-12)
    P_mean_o = Q_o_dim / (π * a_o_dim * b_o_dim + 1e-12)
    
    F_roll_i_dim, F_side_i_dim, H_slide_i_dim, H_spin_i_dim, M_sp_i_dim = integrate_tehd_contact_force(
        u_roll_i_dim, u_side_i_dim, P_mean_i, a_i_dim, b_i_dim, u_entrain_i,
        ω_spin_i, ω_roll_abs_i, R_ball_dim, p.geom.f_i,
        p.trac.A, p.trac.B, p.trac.C, p.trac.D, p.lub.Lambda_LSS, p.lub.beta_temp, p.lub.effusivity
    )
    
    F_roll_o_dim, F_side_o_dim, H_slide_o_dim, H_spin_o_dim, M_sp_o_dim = integrate_tehd_contact_force(
        u_roll_o_dim, u_side_o_dim, P_mean_o, a_o_dim, b_o_dim, u_entrain_o,
        ω_spin_o, ω_roll_abs_o, R_ball_dim, p.geom.f_o,
        p.trac.A, p.trac.B, p.trac.C, p.trac.D, p.lub.Lambda_LSS, p.lub.beta_temp, p.lub.effusivity
    )
    
    F_trac_tang_i = (F_roll_i_dim * ctx.entrain_factor) / ctx.Q_scale
    F_trac_tang_o = (F_roll_o_dim * ctx.entrain_factor) / ctx.Q_scale
    F_side_i = (F_side_i_dim * ctx.entrain_factor) / ctx.Q_scale
    F_side_o = (F_side_o_dim * ctx.entrain_factor) / ctx.Q_scale
    M_spin_i_nd = (M_sp_i_dim * ctx.entrain_factor) / (ctx.Q_scale * ctx.L_scale)
    M_spin_o_nd = (M_sp_o_dim * ctx.entrain_factor) / (ctx.Q_scale * ctx.L_scale)
    
    P_slide_i_nd = (H_slide_i_dim * ctx.entrain_factor) / (ctx.Q_scale * ctx.V_scale)
    P_slide_o_nd = (H_slide_o_dim * ctx.entrain_factor) / (ctx.Q_scale * ctx.V_scale)
    P_spin_i_nd  = (H_spin_i_dim * ctx.entrain_factor) / (ctx.Q_scale * ctx.V_scale)
    P_spin_o_nd  = (H_spin_o_dim * ctx.entrain_factor) / (ctx.Q_scale * ctx.V_scale)
    
    H_i_tot = P_slide_i_nd + P_spin_i_nd
    H_o_tot = P_slide_o_nd + P_spin_o_nd
    
    v_rel_θ_nd = r_b * (θ̇_b - p.load.omega_cage)
    V_ball_dim = sqrt(ẋ_b^2 + v_rel_θ_nd^2 + ṙ_b^2 + 1e-30) * ctx.V_scale
    F_drag_dim = drag_force(V_ball_dim, D_star * ctx.L_scale, p.lub.rho_eff, ctx.μ_oil_dyn, p.cage.cage_web * ctx.L_scale)
    F_drag_total = F_drag_dim / ctx.Q_scale
    F_drag_x = -F_drag_total * (ẋ_b * ctx.V_scale / V_ball_dim)
    F_drag_r = -F_drag_total * (ṙ_b * ctx.V_scale / V_ball_dim)
    F_drag_θ = -F_drag_total * (v_rel_θ_nd * ctx.V_scale / V_ball_dim)
    
    ω_mag_dim = sqrt(ω_body[1]^2 + ω_body[2]^2 + ω_body[3]^2 + 1e-30) * ctx.W_scale
    M_ch_dim = churning_moment(ω_mag_dim, R_ball_k * ctx.L_scale, p.lub.rho_eff, ctx.μ_oil_dyn)
    M_ch_nondim = M_ch_dim / (ctx.Q_scale * ctx.L_scale)
    
    P_drag_nd = F_drag_total * (V_ball_dim / ctx.V_scale)
    P_drag_spin_nd = M_ch_nondim * (ω_mag_dim / ctx.W_scale)
    
    D_ball_dim = D_star * ctx.L_scale
    R_half_D = T_u(0.5) * D_ball_dim
    ndm_fac = T_u(60000.0 / π)
    
    U_ehl_i = p.lub.mu_0 * u_entrain_i / (ctx.E_prime_dim * R_half_D + T_u(1e-30))
    W_ehl_i = Q_i_dim / (ctx.E_prime_dim * R_half_D^2 + T_u(1e-30))
    F_ehl_i = T_u(4.318) * ctx.E_prime_dim * R_half_D^2 / ctx.G_val * (U_ehl_i * ctx.G_val + T_u(1e-12))^T_u(0.658) * (W_ehl_i + T_u(1e-12))^T_u(0.0126)
    φ_sh_i = T_u(1.0) / (T_u(1.0) + T_u(1.84e-9) * (u_entrain_i * ndm_fac + T_u(1e-12))^T_u(1.28) * ctx.ν_cst_064)
    F_roll_ehl_i = F_ehl_i * φ_sh_i
    
    U_ehl_o = p.lub.mu_0 * u_entrain_o / (ctx.E_prime_dim * R_half_D + T_u(1e-30))
    W_ehl_o = Q_o_dim / (ctx.E_prime_dim * R_half_D^2 + T_u(1e-30))
    F_ehl_o = T_u(4.318) * ctx.E_prime_dim * R_half_D^2 / ctx.G_val * (U_ehl_o * ctx.G_val + T_u(1e-12))^T_u(0.658) * (W_ehl_o + T_u(1e-12))^T_u(0.0126)
    φ_sh_o = T_u(1.0) / (T_u(1.0) + T_u(1.84e-9) * (u_entrain_o * ndm_fac + T_u(1e-12))^T_u(1.28) * ctx.ν_cst_064)
    F_roll_ehl_o = F_ehl_o * φ_sh_o
    
    P_ehl_roll_nd = (F_roll_ehl_i * u_entrain_i + F_roll_ehl_o * u_entrain_o) / (ctx.Q_scale * ctx.V_scale)
    
    T_scale_inv = ctx.V_scale / ctx.L_scale
    tau_pass_i = conservative_contact_pass_time(a_i_dim, b_i_dim, u_entrain_i) * T_scale_inv
    tau_pass_o = conservative_contact_pass_time(a_o_dim, b_o_dim, u_entrain_o) * T_scale_inv
    tau_a2_i = tau_pass_i + 0.0
    tau_a2_o = tau_pass_o + 0.0
    
    alpha2_kin_i = asin(smooth_signed_saturate(u_side_i_dim / (u_entrain_i + 1e-6), 0.99))
    alpha2_kin_o = asin(smooth_signed_saturate(u_side_o_dim / (u_entrain_o + 1e-6), 0.99))
    
    return (
        # ── Contact Geometry ──
        dx_i=dx_i, dr_i=dr_i, dx_o=dx_o, dr_o=dr_o,
        δ_i_sm=δ_i_sm, δ_o_sm=δ_o_sm,
        a_i_dim=a_i_dim, b_i_dim=b_i_dim, a_o_dim=a_o_dim, b_o_dim=b_o_dim,
        R_ball_k=R_ball_k,
        X_rel_i=X_rel_i, Y_rel_i=Y_rel_i, Z_rel_i=Z_rel_i,
        n_i_loc=n_i_loc, t_lat_i_loc=t_lat_i_loc, t_roll_i_loc=t_roll_i_loc,
        n_o_loc=n_o_loc, t_lat_o_loc=t_lat_o_loc, t_roll_o_loc=t_roll_o_loc,

        # ── Contact Forces (nondimensional for ODE) ──
        Q_i=Q_i, Q_o=Q_o,
        F_trac_tang_i=F_trac_tang_i, F_trac_tang_o=F_trac_tang_o,
        F_side_i=F_side_i, F_side_o=F_side_o,
        M_spin_i_nd=M_spin_i_nd, M_spin_o_nd=M_spin_o_nd,
        F_drag_x=F_drag_x, F_drag_r=F_drag_r, F_drag_θ=F_drag_θ,
        M_ch_nondim=M_ch_nondim,
        F_roll_ehl_i_nd=F_roll_ehl_i / ctx.Q_scale, F_roll_ehl_o_nd=F_roll_ehl_o / ctx.Q_scale,

        # ── Kinematics (velocities, entrainment) ──
        u_roll_i_dim=u_roll_i_dim, u_roll_o_dim=u_roll_o_dim,
        u_side_i_dim=u_side_i_dim, u_side_o_dim=u_side_o_dim,
        u_entrain_i=u_entrain_i, u_entrain_o=u_entrain_o,
        v_ball_θ=v_ball_θ, v_ir_θ_dim=v_ir_θ_dim, v_or_θ_dim=v_or_θ_dim,

        # ── Power & Thermal (nondimensional for thermal network) ──
        H_i_tot=H_i_tot, H_o_tot=H_o_tot,
        P_drag_nd=P_drag_nd, P_drag_spin_nd=P_drag_spin_nd, P_ehl_roll_nd=P_ehl_roll_nd,

        # ── Dimensional outputs (for field output / postprocessing) ──
        Q_i_dim=Q_i_dim, Q_o_dim=Q_o_dim,
        ω_spin_i_dim=ω_spin_i, ω_spin_o_dim=ω_spin_o,
        F_trac_i_dim=F_roll_i_dim * ctx.entrain_factor, F_trac_o_dim=F_roll_o_dim * ctx.entrain_factor,
        M_sp_i_dim=M_sp_i_dim * ctx.entrain_factor, M_sp_o_dim=M_sp_o_dim * ctx.entrain_factor,
        H_slide_i_dim=H_slide_i_dim * ctx.entrain_factor, H_slide_o_dim=H_slide_o_dim * ctx.entrain_factor,
        H_spin_i_dim=H_spin_i_dim * ctx.entrain_factor, H_spin_o_dim=H_spin_o_dim * ctx.entrain_factor,
        F_drag_dim=F_drag_dim, M_ch_dim=M_ch_dim,
        ω_ball_x_dim=ω_body[1] * ctx.W_scale, ω_ball_y_dim=ω_body[2] * ctx.W_scale, ω_ball_z_dim=ω_body[3] * ctx.W_scale
    )
end
