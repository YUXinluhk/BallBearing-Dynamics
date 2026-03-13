# =====================================================================
# Config/types.jl — All configuration structs for ADORE
#
# Julia port of bearing.py, ehl_traction.py, cage.py, churning.py,
# dynamics.py::SimulationConfig, thermal.py, integrator.py.
#
# Design: immutable `@kwdef` structs + free functions for derived quantities.
# All SI units: [m], [N], [Pa], [rad], [s].
# =====================================================================

using LinearAlgebra

# ── Material ──────────────────────────────────────────────────────────

"""
    MaterialParams(; E, nu)

Elastic material properties for contacting bodies.
"""
Base.@kwdef struct MaterialParams
    E::Float64     # Young's modulus [Pa]
    nu::Float64    # Poisson's ratio [-]
    effusivity::Float64 = 12000.0  # Thermal effusivity [W·s^0.5/(m^2·K)]

    function MaterialParams(E, nu, effusivity=12000.0)
        check_positive("E (Young's modulus)", E, "MaterialParams")
        check_range("nu (Poisson's ratio)", nu, 0.0, 0.5, "MaterialParams";
            inclusive=:neither)
        warn_range("E (Young's modulus)", E, 1e9, 1e12, "MaterialParams")
        check_positive("effusivity", effusivity, "MaterialParams")
        new(E, nu, effusivity)
    end
end

"Plane strain modulus: E/(1 - ν²)"
E_plane_strain(m::MaterialParams) = m.E / (1.0 - m.nu^2)

"Hertzian reduced modulus for identical materials: E/(2(1 - ν²))"
E_star(m::MaterialParams) = m.E / (2.0 * (1.0 - m.nu^2))

"""
    composite_modulus(m1, m2)

Composite elastic modulus E* for two contacting bodies.
1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂
"""
function composite_modulus(m1::MaterialParams, m2::MaterialParams)
    inv_E = (1.0 - m1.nu^2) / m1.E + (1.0 - m2.nu^2) / m2.E
    return 1.0 / inv_E
end

# ── Bearing Geometry ──────────────────────────────────────────────────

"""
    BearingGeometry(; d, n_balls, f_i, f_o, d_m, alpha_0, P_d=0, rho_ball=7800)

Angular contact ball bearing geometry. All lengths [m], angles [rad].
"""
Base.@kwdef struct BearingGeometry
    d::Float64           # Ball diameter [m]
    n_balls::Int         # Number of balls
    f_i::Float64         # Inner groove conformity r_i/d
    f_o::Float64         # Outer groove conformity r_o/d
    d_m::Float64         # Pitch diameter [m]
    alpha_0::Float64     # Free contact angle [rad]
    P_d::Float64 = 0.0   # Diametral clearance [m]
    rho_ball::Float64 = 7800.0  # Ball density [kg/m³]

    function BearingGeometry(d, n_balls, f_i, f_o, d_m, alpha_0, P_d, rho_ball)
        check_positive("d (ball diameter)", d, "BearingGeometry")
        n_balls < 3 && throw(ArgumentError(
            "n_balls (BearingGeometry) must be >= 3, got $n_balls"))
        check_range("f_i (inner conformity)", f_i, 0.5, 0.6, "BearingGeometry";
            inclusive=:neither)
        check_range("f_o (outer conformity)", f_o, 0.5, 0.6, "BearingGeometry";
            inclusive=:neither)
        check_positive("d_m (pitch diameter)", d_m, "BearingGeometry")
        d >= d_m && throw(ArgumentError(
            "d ($d) must be < d_m ($d_m)"))
        d / d_m >= 0.5 && throw(ArgumentError(
            "γ' = d/d_m = $(d/d_m) must be < 0.5"))
        check_range("alpha_0 (contact angle)", alpha_0, 0.0, π / 3, "BearingGeometry")
        check_non_negative("P_d (clearance)", P_d, "BearingGeometry")
        check_positive("rho_ball (density)", rho_ball, "BearingGeometry")
        warn_range("rho_ball (density)", rho_ball, 1000.0, 20000.0, "BearingGeometry")
        new(d, n_balls, f_i, f_o, d_m, alpha_0, P_d, rho_ball)
    end
end

# ── BearingGeometry derived quantities (free functions) ──

r_ball(g::BearingGeometry) = g.d / 2
r_i(g::BearingGeometry) = g.f_i * g.d
r_o(g::BearingGeometry) = g.f_o * g.d

"""
    alpha_free(g::BearingGeometry)

Free contact angle accounting for *additional* diametral clearance P_d.
Harris §7.2 adapted for angular-contact bearings where α₀ is the
designed contact angle at zero additional clearance:

    cos(α_f) = cos(α₀) - P_d / (2·B_D)

where B_D = (f_i - 0.5)·d + (f_o - 0.5)·d.
When P_d = 0, returns α₀.  α_f increases monotonically with P_d.
"""
function alpha_free(g::BearingGeometry)
    if g.P_d <= 0.0
        return g.alpha_0
    end
    B_D = (g.f_i - 0.5) * g.d + (g.f_o - 0.5) * g.d
    cos_af = cos(g.alpha_0) - g.P_d / (2.0 * B_D)
    cos_af = clamp(cos_af, -1.0, 1.0)  # guard
    return acos(cos_af)
end

"Inner race groove center diameter: D_i = d_m + (2f_i - 1)·d·cos(α_f)"
D_i(g::BearingGeometry) = g.d_m + (2g.f_i - 1) * g.d * cos(alpha_free(g))

"Outer race groove center diameter: D_o = d_m - (2f_o - 1)·d·cos(α_f)"
D_o(g::BearingGeometry) = g.d_m - (2g.f_o - 1) * g.d * cos(alpha_free(g))

"Ball mass [kg] — solid sphere"
ball_mass(g::BearingGeometry) = g.rho_ball * (4 / 3) * π * r_ball(g)^3

"Ball moment of inertia [kg·m²] about any diameter: (2/5)mr²"
ball_inertia(g::BearingGeometry) = 0.4 * ball_mass(g) * r_ball(g)^2

"Inner race mass [kg] — annular ring estimate"
function inner_race_mass(g::BearingGeometry)
    R_outer = g.d_m / 2
    R_inner = R_outer - g.d * 0.6
    width = g.d * 1.2
    vol = π * (R_outer^2 - R_inner^2) * width
    return g.rho_ball * vol
end

"Angular spacing between ball centers [rad]"
ball_spacing(g::BearingGeometry) = 2π / g.n_balls

"Σρ̄ for ball/inner-race contact"
function sum_rho_inner(g::BearingGeometry)
    cos_a = cos(alpha_free(g))
    ρ_ball = 2.0 / g.d
    ρ_ir_roll = 2.0 * cos_a / (g.d_m - g.d * cos_a)
    ρ_ir_cross = -1.0 / r_i(g)
    return 2ρ_ball + ρ_ir_roll + ρ_ir_cross
end

"Σρ̄ for ball/outer-race contact"
function sum_rho_outer(g::BearingGeometry)
    cos_a = cos(alpha_free(g))
    ρ_ball = 2.0 / g.d
    ρ_or_roll = -2.0 * cos_a / (g.d_m + g.d * cos_a)
    ρ_or_cross = -1.0 / r_o(g)
    return 2ρ_ball + ρ_or_roll + ρ_or_cross
end

"Curvature difference F(ρ) for ball/inner contact"
function F_rho_inner(g::BearingGeometry)
    cos_a = cos(alpha_free(g))
    ρ_ir_roll = 2.0 * cos_a / (g.d_m - g.d * cos_a)
    ρ_ir_cross = -1.0 / r_i(g)
    sr = sum_rho_inner(g)
    abs(sr) < 1e-30 && return 0.0
    return abs(ρ_ir_cross - ρ_ir_roll) / sr
end

"Curvature difference F(ρ) for ball/outer contact"
function F_rho_outer(g::BearingGeometry)
    cos_a = cos(alpha_free(g))
    ρ_or_roll = -2.0 * cos_a / (g.d_m + g.d * cos_a)
    ρ_or_cross = -1.0 / r_o(g)
    sr = sum_rho_outer(g)
    abs(sr) < 1e-30 && return 0.0
    return abs(ρ_or_cross - ρ_or_roll) / sr
end

# ── Standard catalog ─────────────────────────────────────────────────

"7210B angular contact ball bearing (standard catalog)"
function bearing_7210B()
    geom = BearingGeometry(
        d=12.7e-3,
        n_balls=16,
        f_i=0.52,
        f_o=0.52,
        d_m=70.0e-3,
        alpha_0=deg2rad(40.0),
        P_d=0.0,
        rho_ball=7800.0,
    )
    mat = MaterialParams(E=2.08e11, nu=0.3)
    return geom, mat
end

"7008C angular contact ball bearing (40mm bore, 15° contact angle)"
function bearing_7008C()
    # ISO 7008C: d_bore=40mm, d_outer=68mm, width=15mm
    # Ball d=7.938mm, Z=16, α₀=15°, f_i=0.52, f_o=0.53
    geom = BearingGeometry(
        d=7.938e-3,
        n_balls=16,
        f_i=0.52,
        f_o=0.53,
        d_m=54.0e-3,   # (40+68)/2
        alpha_0=deg2rad(15.0),
        P_d=0.0,
        rho_ball=7800.0,
    )
    mat = MaterialParams(E=2.08e11, nu=0.3)
    return geom, mat
end

"7010C angular contact ball bearing (50mm bore, 15° contact angle)"
function bearing_7010C()
    # ISO 7010C: d_bore=50mm, d_outer=80mm, width=16mm
    # Ball d=9.525mm, Z=16, α₀=15°, f_i=0.52, f_o=0.53
    geom = BearingGeometry(
        d=9.525e-3,
        n_balls=16,
        f_i=0.52,
        f_o=0.53,
        d_m=65.0e-3,   # (50+80)/2
        alpha_0=deg2rad(15.0),
        P_d=0.0,
        rho_ball=7800.0,
    )
    mat = MaterialParams(E=2.08e11, nu=0.3)
    return geom, mat
end

"7014C angular contact ball bearing (70mm bore, 15° contact angle)"
function bearing_7014C()
    # ISO 7014C: d_bore=70mm, d_outer=110mm, width=20mm
    # Ball d=12.7mm, Z=18, α₀=15°, f_i=0.52, f_o=0.53
    geom = BearingGeometry(
        d=12.7e-3,
        n_balls=18,
        f_i=0.52,
        f_o=0.53,
        d_m=90.0e-3,   # (70+110)/2
        alpha_0=deg2rad(15.0),
        P_d=0.0,
        rho_ball=7800.0,
    )
    mat = MaterialParams(E=2.08e11, nu=0.3)
    return geom, mat
end

"""
    bearing_custom(; bore, outer, d_ball, n_balls, alpha_deg, f_i=0.52, f_o=0.53)

Create a custom bearing from basic catalog dimensions.
"""
function bearing_custom(; bore, outer, d_ball, n_balls, alpha_deg, f_i=0.52, f_o=0.53)
    geom = BearingGeometry(
        d=d_ball,
        n_balls=n_balls,
        f_i=f_i,
        f_o=f_o,
        d_m=(bore + outer) / 2.0,
        alpha_0=deg2rad(alpha_deg),
        P_d=0.0,
        rho_ball=7800.0,
    )
    mat = MaterialParams(E=2.08e11, nu=0.3)
    return geom, mat
end

# ── Lubricant ─────────────────────────────────────────────────────────

"""
    LubricantParams(; mu_0, alpha_pv, beta_temp=0.05, ...)

EHL lubricant properties. All SI units.
"""
Base.@kwdef struct LubricantParams
    mu_0::Float64           # Reference dynamic viscosity [Pa·s]
    alpha_pv::Float64       # Pressure-viscosity coefficient [1/Pa]
    beta_temp::Float64 = 0.05   # Temperature-viscosity coefficient [1/K]
    Lambda_LSS::Float64 = 0.04  # Limiting shear stress coefficient [-]
    T_0::Float64 = 350.0       # Reference temperature [K]
    K_th::Float64 = 0.14       # Thermal conductivity [W/(m·K)]
    rho_lub::Float64 = 860.0   # Lubricant density [kg/m³]
    c_p::Float64 = 2000.0      # Specific heat capacity [J/(kg·K)]

    function LubricantParams(mu_0, alpha_pv, beta_temp, Lambda_LSS, T_0, K_th,
        rho_lub, c_p)
        check_positive("mu_0 (viscosity)", mu_0, "LubricantParams")
        check_positive("alpha_pv", alpha_pv, "LubricantParams")
        check_range("T_0 (reference temp)", T_0, 250.0, 500.0, "LubricantParams")
        check_positive("K_th", K_th, "LubricantParams")
        check_positive("rho_lub", rho_lub, "LubricantParams")
        new(mu_0, alpha_pv, beta_temp, Lambda_LSS, T_0, K_th, rho_lub, c_p)
    end
end

"Standard MIL-L-23699 lubricant"
function lubricant_mil_l_23699()
    LubricantParams(
        mu_0=0.005,
        alpha_pv=1.5e-8,
        beta_temp=0.04,
        Lambda_LSS=0.04,
        T_0=350.0,
        K_th=0.14,
    )
end

# ── Traction ──────────────────────────────────────────────────────────

"""
    TractionParams(; kappa_0, kappa_inf, kappa_m, u_m)

4-parameter traction curve: κ = (A + B·u)·exp(-C·u) + D
Boundary conditions solved at construction time (A, B, C, D stored).
"""
struct TractionParams
    kappa_0::Float64    # Initial traction coefficient κ(0)
    kappa_inf::Float64  # Asymptotic traction coefficient κ(∞)
    kappa_m::Float64    # Peak traction coefficient
    u_m::Float64        # Slip velocity at peak [m/s]
    A::Float64          # Solved parameter
    B::Float64          # Solved parameter
    C::Float64          # Solved parameter
    D::Float64          # Solved parameter

    function TractionParams(; kappa_0, kappa_inf, kappa_m, u_m)
        check_positive("u_m", u_m, "TractionParams")
        kappa_inf >= kappa_0 && throw(ArgumentError(
            "kappa_inf ($kappa_inf) must be < kappa_0 ($kappa_0)"))

        # Solve for B from transcendental equation (Eq 4.8)
        # (κ₀ - κ∞ + B·u_m)·exp(-B·u_m/(κ₀ - κ∞ + B·u_m)) = κ_m - κ∞
        Δκ = kappa_0 - kappa_inf
        target = kappa_m - kappa_inf

        function residual(B_val)
            x = Δκ + B_val * u_m
            x <= 0 && return -target  # guard
            return x * exp(-B_val * u_m / x) - target
        end

        # Bisection search for B (matching Python: B_max = 5 * target / u_m)
        B_lo, B_hi = 0.0, 5.0 * target / u_m
        for _ in 1:100
            B_mid = 0.5 * (B_lo + B_hi)
            if residual(B_mid) > 0
                B_hi = B_mid
            else
                B_lo = B_mid
            end
        end
        B_val = 0.5 * (B_lo + B_hi)

        # Derive A, C, D from B
        D = kappa_inf
        A = kappa_0 - D
        C = B_val / (A + B_val * u_m)  # from peak condition: dκ/du = 0 at u_m

        new(kappa_0, kappa_inf, kappa_m, u_m, A, B_val, C, D)
    end
end

"Default traction parameters"
function traction_params_default()
    TractionParams(kappa_0=0.08, kappa_inf=0.04, kappa_m=0.10, u_m=0.5)
end

# ── Cage Geometry ─────────────────────────────────────────────────────

"""
    CageGeometry(; pocket_radius, pocket_clearance, n_pockets, ...)

Cage geometry parameters for ball/cage interaction.
"""
Base.@kwdef struct CageGeometry
    pocket_radius::Float64       # [m]
    pocket_clearance::Float64    # [m]
    n_pockets::Int
    cage_inner_radius::Float64   # [m]
    cage_outer_radius::Float64   # [m]
    pilot_clearance::Float64     # [m]
    pilot_is_inner::Bool = true
    cage_mass::Float64           # [kg]
    cage_inertia_xx::Float64     # [kg·m²]
    cage_inertia_yy::Float64     # [kg·m²]
    stiffness_pocket::Float64 = 1e8     # [N/m]
    stiffness_pilot::Float64 = 1e9     # [N/m]
    mu_pocket::Float64 = 0.05
    mu_pilot::Float64 = 0.03
    c_damping::Float64 = 50.0   # [N·s/m]

    function CageGeometry(pocket_radius, pocket_clearance, n_pockets,
        cage_inner_radius, cage_outer_radius, pilot_clearance,
        pilot_is_inner, cage_mass, cage_inertia_xx, cage_inertia_yy,
        stiffness_pocket, stiffness_pilot,
        mu_pocket, mu_pilot, c_damping)
        check_positive("pocket_radius", pocket_radius, "CageGeometry")
        check_positive("pocket_clearance", pocket_clearance, "CageGeometry")
        n_pockets < 3 && throw(ArgumentError("n_pockets must be >= 3"))
        check_positive("cage_mass", cage_mass, "CageGeometry")
        new(pocket_radius, pocket_clearance, n_pockets,
            cage_inner_radius, cage_outer_radius, pilot_clearance,
            pilot_is_inner, cage_mass, cage_inertia_xx, cage_inertia_yy,
            stiffness_pocket, stiffness_pilot, mu_pocket, mu_pilot, c_damping)
    end
end

"Create CageGeometry from BearingGeometry with default ratios"
function cage_from_bearing(g::BearingGeometry;
    pocket_clearance_ratio=0.01,
    pilot_clearance_ratio=0.002)
    rb = r_ball(g)
    pocket_r = rb * 1.01   # slightly larger than ball
    pocket_clr = g.d * pocket_clearance_ratio

    # Cage radii: simple approximation
    cage_ir = g.d_m / 2 - rb * 1.2
    cage_or = g.d_m / 2 + rb * 1.2
    pilot_clr = g.d_m * pilot_clearance_ratio

    # Cage mass: thin annular ring
    cage_width = g.d * 1.0
    cage_thickness = rb * 0.4
    R_mid = g.d_m / 2
    cage_vol = 2π * R_mid * cage_thickness * cage_width
    ρ_cage = 2700.0   # aluminium typical
    cage_m = ρ_cage * cage_vol

    # Cage inertia (thin ring approximation)
    I_xx = cage_m * R_mid^2
    I_yy = 0.5 * cage_m * R_mid^2

    CageGeometry(
        pocket_radius=pocket_r,
        pocket_clearance=pocket_clr,
        n_pockets=g.n_balls,
        cage_inner_radius=cage_ir,
        cage_outer_radius=cage_or,
        pilot_clearance=pilot_clr,
        pilot_is_inner=true,
        cage_mass=cage_m,
        cage_inertia_xx=I_xx,
        cage_inertia_yy=I_yy,
    )
end

# ── Churning Parameters ──────────────────────────────────────────────

"""
    ChurningParams(; rho_oil=860, rho_air=0.6, mu_oil=5e-3, fill_fraction=0.02, cage_web=0.5e-3)
"""
Base.@kwdef struct ChurningParams
    rho_oil::Float64 = 860.0
    rho_air::Float64 = 0.6
    mu_oil::Float64 = 5.0e-3
    fill_fraction::Float64 = 0.02
    cage_web::Float64 = 0.5e-3

    function ChurningParams(rho_oil, rho_air, mu_oil, fill_fraction, cage_web)
        check_positive("rho_oil", rho_oil, "ChurningParams")
        check_non_negative("rho_air", rho_air, "ChurningParams")
        check_positive("mu_oil", mu_oil, "ChurningParams")
        check_range("fill_fraction", fill_fraction, 0.0, 1.0, "ChurningParams")
        check_non_negative("cage_web", cage_web, "ChurningParams")
        new(rho_oil, rho_air, mu_oil, fill_fraction, cage_web)
    end
end

"Effective density of oil-air mixture"
rho_eff(p::ChurningParams) = p.rho_oil * p.fill_fraction + p.rho_air * (1 - p.fill_fraction)

# ── Integrator Config ─────────────────────────────────────────────────

"""
    IntegratorConfig(; rtol=1e-6, atol=1e-9, ...)
"""
Base.@kwdef struct IntegratorConfig
    rtol::Float64 = 1e-6
    atol::Float64 = 1e-9
    h_max::Float64 = 1e-3
    max_steps::Int = 100_000
end

# ── Thermal Parameters ────────────────────────────────────────────────

"""
    ThermalParams(; T_ambient=300, ...)

Lumped Parameter Thermal Network (LPTN) 4-node model parameters.
"""
Base.@kwdef struct ThermalParams
    T_ambient::Float64 = 300.0       # [K]
    T_init::Float64 = 350.0          # [K]
    C_ir::Float64 = NaN              # inner race heat capacity [J/K], NaN → auto
    C_or::Float64 = NaN              # outer race
    C_ball::Float64 = NaN            # all balls
    C_oil::Float64 = NaN             # oil sump
    G_ir_ball::Float64 = 50.0        # [W/K]
    G_or_ball::Float64 = 50.0
    G_ball_oil::Float64 = 20.0
    G_or_amb::Float64 = 10.0
    G_oil_amb::Float64 = 5.0
    CTE::Float64 = 12.0e-6          # [1/K]
    T_ref::Float64 = 293.15         # [K]
    eta_ball::Float64 = 0.5
    c_steel::Float64 = 460.0        # [J/(kg·K)]
    oil_flow_rate::Float64 = 0.0    # [cm³/min] oil flow rate (0 = no circulation)
    T_oil_inlet::Float64 = NaN      # [K] oil inlet temperature (NaN → use T_init)
end

"Mutable LPTN 4-node temperature state"
mutable struct LPTNState
    T_ir::Float64
    T_or::Float64
    T_ball::Float64
    T_oil::Float64
end

LPTNState() = LPTNState(350.0, 350.0, 350.0, 350.0)
LPTNState(T::Float64) = LPTNState(T, T, T, T)

as_array(s::LPTNState) = @SVector [s.T_ir, s.T_or, s.T_ball, s.T_oil]

function lptn_from_array(arr)
    LPTNState(arr[1], arr[2], arr[3], arr[4])
end

Base.copy(s::LPTNState) = LPTNState(s.T_ir, s.T_or, s.T_ball, s.T_oil)

# ── Simulation Config ─────────────────────────────────────────────────

"""
    SimulationConfig(; t_end, dt_output, inner_race_speed, ...)
"""
Base.@kwdef struct SimulationConfig
    t_end::Float64 = 0.01
    dt_output::Float64 = 1e-5
    inner_race_speed::Float64 = 0.0     # [rad/s]
    outer_race_speed::Float64 = 0.0     # [rad/s]
    F_axial::Float64 = 0.0              # [N]
    F_radial::Float64 = 0.0             # [N]
    t_ramp_end::Float64 = 0.0           # [s]
    mu_spin::Float64 = 0.06             # spin friction coefficient
    c_structural::Float64 = 10.0        # structural damping [N·s/m]
    zeta::Float64 = 0.03                # damping ratio
    delta_r_thermal::Float64 = 0.0      # thermal radial expansion [m]
    integrator::IntegratorConfig = IntegratorConfig()
    churning::ChurningParams = ChurningParams()
    thermal::ThermalParams = ThermalParams()

    function SimulationConfig(t_end, dt_output, inner_race_speed, outer_race_speed, F_axial, F_radial,
        t_ramp_end, mu_spin, c_structural, zeta,
        delta_r_thermal, integrator, churning, thermal)
        check_positive("t_end", t_end, "SimulationConfig")
        dt_output >= t_end && throw(ArgumentError(
            "dt_output ($dt_output) must be < t_end ($t_end)"))
        check_non_negative("inner_race_speed", inner_race_speed, "SimulationConfig")
        check_range("zeta", zeta, 0.0, 1.0, "SimulationConfig")
        new(t_end, dt_output, inner_race_speed, outer_race_speed, F_axial, F_radial,
            t_ramp_end, mu_spin, c_structural, zeta, delta_r_thermal,
            integrator, churning, thermal)
    end
end
