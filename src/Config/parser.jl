# =====================================================================
# Config/parser.jl — TOML Configuration Parser for ADORE
# =====================================================================

using TOML

"""
    load_simulation_config(filepath::String)

Reads a TOML configuration file and constructs the 6 core simulation parameter structs:
`geom`, `mat`, `lub`, `trac`, `cage`, `config`
"""
function load_simulation_config(filepath::String)
    data = TOML.parsefile(filepath)

    # ── 1. Geometry ──
    if haskey(data, "geometry") && haskey(data["geometry"], "type")
        # Load from external TOML catalog if type is given
        gtype = string(data["geometry"]["type"])
        catalog_path = normpath(joinpath(@__DIR__, "../../data/catalogs/bearings.toml"))
        if !isfile(catalog_path)
            error("Bearing catalog not found at $catalog_path")
        end
        catalogs = TOML.parsefile(catalog_path)
        if !haskey(catalogs, gtype)
            error("Unknown bearing geometry preset '$gtype' not found in catalog")
        end
        bg_cat = catalogs[gtype]["geometry"]
        geom = BearingGeometry(
            d=bg_cat["d"], n_balls=bg_cat["n_balls"],
            f_i=bg_cat["f_i"], f_o=bg_cat["f_o"],
            d_m=bg_cat["d_m"], alpha_0=bg_cat["alpha_0"],
            P_d=get(bg_cat, "P_d", 0.0), rho_ball=get(bg_cat, "rho_ball", 7800.0)
        )
        if haskey(catalogs[gtype], "material")
            mat_cat = catalogs[gtype]["material"]
            mat = MaterialParams(E=mat_cat["E"], nu=mat_cat["nu"])
        else
            mat = MaterialParams(E=2.08e11, nu=0.3)
        end
        # Overwrite specific keys if they are present
        bg_args = Dict{Symbol,Any}()
        for key in fieldnames(BearingGeometry)
            skey = string(key)
            if haskey(data["geometry"], skey)
                bg_args[key] = data["geometry"][skey]
            else
                bg_args[key] = getproperty(geom, key)
            end
        end
        geom = BearingGeometry(bg_args[:d], bg_args[:n_balls], bg_args[:f_i], bg_args[:f_o], bg_args[:d_m], bg_args[:alpha_0], bg_args[:P_d], bg_args[:rho_ball])
    else
        bg = data["geometry"]
        geom = BearingGeometry(
            d=bg["d"],
            n_balls=bg["n_balls"],
            f_i=bg["f_i"],
            f_o=bg["f_o"],
            d_m=bg["d_m"],
            alpha_0=bg["alpha_0"],
            P_d=get(bg, "P_d", 0.0),
            rho_ball=get(bg, "rho_ball", 7800.0)
        )
        mat = MaterialParams(E=2.08e11, nu=0.3)
    end

    # ── 2. Material ──
    if haskey(data, "material")
        dm = data["material"]
        mat = MaterialParams(
            E=get(dm, "E", mat.E),
            nu=get(dm, "nu", mat.nu)
        )
    end

    # ── 3. Lubricant ──
    if haskey(data, "lubricant") && get(data["lubricant"], "type", "") == "mil_l_23699"
        lub = lubricant_mil_l_23699()
    else
        lub_data = get(data, "lubricant", Dict())
        # Default to MIL-L-23699 baseline if not fully specified
        default_lub = lubricant_mil_l_23699()
        lub = LubricantParams(
            mu_0=get(lub_data, "mu_0", default_lub.mu_0),
            alpha_pv=get(lub_data, "alpha_pv", default_lub.alpha_pv),
            beta_temp=get(lub_data, "beta_temp", default_lub.beta_temp),
            Lambda_LSS=get(lub_data, "Lambda_LSS", default_lub.Lambda_LSS),
            T_0=get(lub_data, "T_0", default_lub.T_0),
            K_th=get(lub_data, "K_th", default_lub.K_th),
            rho_lub=get(lub_data, "rho_lub", default_lub.rho_lub),
            c_p=get(lub_data, "c_p", default_lub.c_p)
        )
    end

    # ── 4. Traction ──
    trac_data = get(data, "traction", Dict())
    default_trac = traction_params_default()
    trac = TractionParams(
        kappa_0=get(trac_data, "kappa_0", default_trac.kappa_0),
        kappa_inf=get(trac_data, "kappa_inf", default_trac.kappa_inf),
        kappa_m=get(trac_data, "kappa_m", default_trac.kappa_m),
        u_m=get(trac_data, "u_m", default_trac.u_m)
    )

    # ── 5. Cage ──
    cage_data = get(data, "cage", Dict())
    if get(cage_data, "auto_generate", true)
        cage = cage_from_bearing(geom;
            pocket_clearance_ratio=get(cage_data, "pocket_clearance_ratio", 0.01),
            pilot_clearance_ratio=get(cage_data, "pilot_clearance_ratio", 0.002))
    else
        cage = CageGeometry(
            pocket_radius=cage_data["pocket_radius"],
            pocket_clearance=cage_data["pocket_clearance"],
            n_pockets=cage_data["n_pockets"],
            cage_inner_radius=cage_data["cage_inner_radius"],
            cage_outer_radius=cage_data["cage_outer_radius"],
            pilot_clearance=cage_data["pilot_clearance"],
            pilot_is_inner=get(cage_data, "pilot_is_inner", true),
            cage_mass=cage_data["cage_mass"],
            cage_inertia_xx=cage_data["cage_inertia_xx"],
            cage_inertia_yy=cage_data["cage_inertia_yy"],
            stiffness_pocket=get(cage_data, "stiffness_pocket", 1e8),
            stiffness_pilot=get(cage_data, "stiffness_pilot", 1e9),
            mu_pocket=get(cage_data, "mu_pocket", 0.05),
            mu_pilot=get(cage_data, "mu_pilot", 0.03),
            c_damping=get(cage_data, "c_damping", 50.0)
        )
    end

    # ── 6. Simulation Config ──
    cfg = get(data, "operating_conditions", Dict())

    # Integrator
    int_data = get(data, "integrator", Dict())
    integrator = IntegratorConfig(
        rtol=get(int_data, "rtol", 1e-4),
        atol=get(int_data, "atol", 1e-7),
        h_max=get(int_data, "h_max", 1e-4),
        max_steps=get(int_data, "max_steps", 10_000_000),
        eps_contact=get(int_data, "eps_contact", 1e-6)
    )

    # Dynamics/Operation
    inner_race_speed = get(cfg, "inner_race_speed_RPM", 0.0) * π / 30
    outer_race_speed = get(cfg, "outer_race_speed_RPM", 0.0) * π / 30

    churn_data = get(data, "churning", Dict())
    churn = ChurningParams(
        rho_oil=get(churn_data, "rho_oil", 860.0),
        rho_air=get(churn_data, "rho_air", 0.6),
        mu_oil=get(churn_data, "mu_oil", 5.0e-3),
        fill_fraction=get(churn_data, "fill_fraction", 0.02),
        cage_web=get(churn_data, "cage_web", 0.5e-3)
    )

    therm_data = get(data, "thermal", Dict())
    thermal = ThermalParams(
        T_ambient=get(therm_data, "T_ambient", 300.0),
        T_init=get(therm_data, "T_init", 350.0),
        T_ref=get(therm_data, "T_ref", 300.0),
        C_ir=get(therm_data, "C_ir", NaN),
        C_or=get(therm_data, "C_or", NaN),
        C_ball=get(therm_data, "C_ball", NaN),
        C_oil=get(therm_data, "C_oil", NaN),
        G_ir_ball=get(therm_data, "G_ir_ball", 50.0),
        G_or_ball=get(therm_data, "G_or_ball", 50.0),
        G_ball_oil=get(therm_data, "G_ball_oil", 20.0),
        G_or_amb=get(therm_data, "G_or_amb", 10.0),
        G_oil_amb=get(therm_data, "G_oil_amb", 5.0),
        oil_flow_rate=get(therm_data, "oil_flow_rate", 0.0),
        T_oil_inlet=get(therm_data, "T_oil_inlet", NaN),
        th_accel=get(therm_data, "th_accel", 10000.0),
    )

    config = SimulationConfig(
        t_end=get(cfg, "t_end", 0.5),
        dt_output=get(cfg, "dt_output", 5e-3),
        inner_race_speed=inner_race_speed,
        outer_race_speed=outer_race_speed,
        F_axial=get(cfg, "F_axial", 0.0),
        F_radial=get(cfg, "F_radial", 0.0),
        t_ramp_end=get(cfg, "t_ramp_end", 0.01),
        mu_spin=get(cfg, "mu_spin", 0.06),
        c_structural=get(cfg, "c_structural", 0.0),
        zeta=get(cfg, "zeta", 0.10),
        alpha2_relax_time=get(cfg, "alpha2_relax_time", 2e-6),
        delta_r_thermal=get(cfg, "delta_r_thermal", 0.0),
        integrator=integrator,
        churning=churn,
        thermal=thermal
    )

    return geom, mat, lub, trac, cage, config
end
