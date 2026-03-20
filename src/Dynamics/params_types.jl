# =====================================================================
# Dynamics/params_types.jl — Type definitions for structured ODE parameters
#
# Separated from params.jl to allow early inclusion before Physics layer.
# =====================================================================

"""
    GeomParams — Nondimensional bearing geometry for ODE kernel.
"""
struct GeomParams
    D::Float64       # ball diameter ★
    d_m::Float64     # pitch diameter ★
    f_i::Float64     # inner conformity
    f_o::Float64     # outer conformity
    alpha_0::Float64 # free contact angle [rad]
    Z::Int           # number of balls
    D_i::Float64     # inner groove center diameter ★
    D_o::Float64     # outer groove center diameter ★
    drb_i::Float64   # distance race-ball center inner ★
    drb_o::Float64   # distance race-ball center outer ★
    x_gi0::Float64   # inner groove center axial offset ★
    x_go0::Float64   # outer groove center axial offset ★
end

"""
    MassParams — Nondimensional masses and inertias.
"""
struct MassParams
    m_ball::Float64
    J_ball::Float64
    m_ir::Float64
    cage_mass::Float64
    cage_Ixx::Float64
    cage_Iyy::Float64
end

"""
    HertzParams — Nondimensional Hertz contact precomputations.
"""
struct HertzParams
    E_prime::Float64
    Y_i::Float64
    Y_o::Float64
    a_i::Float64
    b_i::Float64
    a_o::Float64
    b_o::Float64
    sr_i::Float64
    sr_o::Float64
    E2_i::Float64
    E2_o::Float64
    eps_contact::Float64
end

"""
    LubParams — Lubricant properties for EHL/traction.
"""
struct LubParams
    mu_0::Float64
    alpha_pv::Float64
    beta_temp::Float64
    Lambda_LSS::Float64
    T_0::Float64
    K_th::Float64
    rho_lub::Float64
    rho_eff::Float64
    mu_oil::Float64
    effusivity::Float64
end

"""
    TractionCoeffs — 4-parameter traction curve coefficients.
"""
struct TractionCoeffs
    A::Float64
    B::Float64
    C::Float64
    D::Float64
end

"""
    DampingCoeffs — Per-mode damping coefficients (nondimensional).
"""
struct DampingCoeffs
    c_ball_trans::Float64
    c_ball_orbit::Float64
    c_ball_spin::Float64
    c_ir::Float64
    c_cage::Float64
    c_tilt::Float64
end

"""
    CageParams — Nondimensional cage geometry and contact.
"""
struct CageParams
    pocket_clr::Float64
    k_pocket::Float64
    mu_pocket::Float64
    k_pilot::Float64
    mu_pilot::Float64
    pilot_clr::Float64
    pilot_is_inner::Bool
    cage_ir::Float64
    cage_or::Float64
    cage_web::Float64
    c_cage::Float64
end

"""
    LoadParams — Applied loads and speeds (nondimensional).
"""
struct LoadParams
    omega_ir::Float64
    omega_or::Float64
    omega_cage::Float64
    F_a::Float64
    F_r::Float64
    t_ramp::Float64
    mu_spin::Float64
    zeta::Float64
end

"""
    ScaleFactors — Dimensional scale factors for kernel use.
"""
struct ScaleFactors
    V::Float64
    W::Float64
    L::Float64
    Q::Float64
    T::Float64
end

"""
    ThermalNDParams — Nondimensional thermal network parameters.
"""
struct ThermalNDParams
    mcp_i::Float64
    mcp_o::Float64
    mcp_ball::Float64
    mcp_oil::Float64
    G_ir_ball::Float64
    G_or_ball::Float64
    G_ball_oil::Float64
    G_or_amb::Float64
    G_oil_amb::Float64
    T_amb::Float64
    CTE::Float64
    T_ref::Float64
    rho_race::Float64
    E_young::Float64
    oil_flow_mdot_cp::Float64
    T_oil_inlet::Float64
end

"""
    ODEParams — Complete structured parameter pack for the ODE kernel.
"""
struct ODEParams
    geom::GeomParams
    mass::MassParams
    hertz::HertzParams
    lub::LubParams
    trac::TractionCoeffs
    damp::DampingCoeffs
    cage::CageParams
    load::LoadParams
    scale::ScaleFactors
    thermal::ThermalNDParams
end
