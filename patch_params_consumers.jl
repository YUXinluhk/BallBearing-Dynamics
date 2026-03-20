#!/usr/bin/env julia
# patch_params_consumers.jl — Replace all p[P_XXX] references with p.sub.field
# This is a comprehensive find-and-replace across 3 major consumer files.

# ── Mapping from old p[P_XXX] to new p.sub.field syntax ──
replacements = [
    # Geometry
    "p[P_D]" => "p.geom.D",
    "p[P_DM]" => "p.geom.d_m",
    "p[P_FI]" => "p.geom.f_i",
    "p[P_FO]" => "p.geom.f_o",
    "p[P_ALPHA0]" => "p.geom.alpha_0",
    "p[P_NBALL]" => "p.geom.Z",
    "p[P_DI]" => "p.geom.D_i",
    "p[P_DO]" => "p.geom.D_o",
    "p[P_DRB_I]" => "p.geom.drb_i",
    "p[P_DRB_O]" => "p.geom.drb_o",
    "p[P_X_GI0]" => "p.geom.x_gi0",
    "p[P_X_GO0]" => "p.geom.x_go0",

    # Mass
    "p[P_MBALL]" => "p.mass.m_ball",
    "p[P_JBALL]" => "p.mass.J_ball",
    "p[P_MIR]" => "p.mass.m_ir",
    "p[P_CAGE_MASS]" => "p.mass.cage_mass",
    "p[P_CAGE_IXX]" => "p.mass.cage_Ixx",
    "p[P_CAGE_IYY]" => "p.mass.cage_Iyy",

    # Hertz
    "p[P_E_PRIME]" => "p.hertz.E_prime",
    "p[P_YI]" => "p.hertz.Y_i",
    "p[P_YO]" => "p.hertz.Y_o",
    "p[P_AI]" => "p.hertz.a_i",
    "p[P_BI]" => "p.hertz.b_i",
    "p[P_AO]" => "p.hertz.a_o",
    "p[P_BO]" => "p.hertz.b_o",
    "p[P_SRI]" => "p.hertz.sr_i",
    "p[P_SRO]" => "p.hertz.sr_o",
    "p[P_E2_I]" => "p.hertz.E2_i",
    "p[P_E2_O]" => "p.hertz.E2_o",
    "p[P_EPS_CONTACT]" => "p.hertz.eps_contact",

    # Lubricant
    "p[P_MU0]" => "p.lub.mu_0",
    "p[P_ALPHA_PV]" => "p.lub.alpha_pv",
    "p[P_BETA_TEMP]" => "p.lub.beta_temp",
    "p[P_LAMBDA_LSS]" => "p.lub.Lambda_LSS",
    "p[P_T0]" => "p.lub.T_0",
    "p[P_KTH]" => "p.lub.K_th",
    "p[P_RHO_LUB]" => "p.lub.rho_lub",
    "p[P_RHO_EFF]" => "p.lub.rho_eff",
    "p[P_MU_OIL]" => "p.lub.mu_oil",
    "p[P_EFFUSIVITY]" => "p.lub.effusivity",

    # Traction
    "p[P_TRAC_A]" => "p.trac.A",
    "p[P_TRAC_B]" => "p.trac.B",
    "p[P_TRAC_C]" => "p.trac.C",
    "p[P_TRAC_D]" => "p.trac.D",

    # Damping
    "p[P_C_BALL_TRANS]" => "p.damp.c_ball_trans",
    "p[P_C_BALL_ORBIT]" => "p.damp.c_ball_orbit",
    "p[P_C_BALL_SPIN]" => "p.damp.c_ball_spin",
    "p[P_C_IR_DAMP]" => "p.damp.c_ir",
    "p[P_C_CAGE_DAMP]" => "p.damp.c_cage",
    "p[P_C_TILT]" => "p.damp.c_tilt",
    "p[P_C_STRUCT]" => "p.damp.c_ir",  # legacy alias
    "p[P_ZETA]" => "p.load.zeta",
    "p[P_MU_SPIN]" => "p.load.mu_spin",

    # Cage
    "p[P_POCKET_CLR]" => "p.cage.pocket_clr",
    "p[P_K_POCKET]" => "p.cage.k_pocket",
    "p[P_MU_POCKET]" => "p.cage.mu_pocket",
    "p[P_K_PILOT]" => "p.cage.k_pilot",
    "p[P_MU_PILOT]" => "p.cage.mu_pilot",
    "p[P_PILOT_CLR]" => "p.cage.pilot_clr",
    "p[P_PILOT_IR]" => "(p.cage.pilot_is_inner ? 1.0 : 0.0)",
    "p[P_CAGE_IR]" => "p.cage.cage_ir",
    "p[P_CAGE_OR]" => "p.cage.cage_or",
    "p[P_CAGE_WEB]" => "p.cage.cage_web",
    "p[P_C_CAGE]" => "p.cage.c_cage",
    "p[P_POCKET_R]" => "p.cage.pocket_clr",  # pocket_radius not used in kernel

    # Loads
    "p[P_OMEGA_IR]" => "p.load.omega_ir",
    "p[P_OMEGA_OR]" => "p.load.omega_or",
    "p[P_OMEGA_CAGE]" => "p.load.omega_cage",
    "p[P_FA]" => "p.load.F_a",
    "p[P_FR]" => "p.load.F_r",
    "p[P_T_RAMP]" => "p.load.t_ramp",

    # Scales
    "p[P_V_SCALE]" => "p.scale.V",
    "p[P_W_SCALE]" => "p.scale.W",
    "p[P_L_SCALE]" => "p.scale.L",
    "p[P_Q_SCALE]" => "p.scale.Q",

    # Thermal
    "p[P_MCP_I]" => "p.thermal.mcp_i",
    "p[P_MCP_O]" => "p.thermal.mcp_o",
    "p[P_MCP_BALL]" => "p.thermal.mcp_ball",
    "p[P_MCP_OIL]" => "p.thermal.mcp_oil",
    "p[P_G_IR_BALL]" => "p.thermal.G_ir_ball",
    "p[P_G_OR_BALL]" => "p.thermal.G_or_ball",
    "p[P_G_BALL_OIL]" => "p.thermal.G_ball_oil",
    "p[P_G_OR_AMB]" => "p.thermal.G_or_amb",
    "p[P_G_OIL_AMB]" => "p.thermal.G_oil_amb",
    "p[P_T_AMB]" => "p.thermal.T_amb",
    "p[P_CTE]" => "p.thermal.CTE",
    "p[P_T_REF]" => "p.thermal.T_ref",
    "p[P_RHO_RACE]" => "p.thermal.rho_race",
    "p[P_E_YOUNG]" => "p.thermal.E_young",
    "p[P_OIL_FLOW_MDOT_CP]" => "p.thermal.oil_flow_mdot_cp",
    "p[P_T_OIL_INLET]" => "p.thermal.T_oil_inlet",
    "p[P_DELTA_R_TH]" => "0.0",  # legacy, unused
    "p[P_ALPHA2_TAU]" => "0.0",  # eliminated
]

function patch_file!(path)
    text = read(path, String)
    changed = false
    for (old, new) in replacements
        if occursin(old, text)
            text = replace(text, old => new)
            changed = true
        end
    end
    
    # Fix the function signature: p::Vector{Float64} → p::ODEParams
    text = replace(text, "p::Vector{Float64}" => "p::ODEParams")
    
    # Fix ctx field: the KinematicContext stores p::Vector{Float64} → p::ODEParams
    # Also fix Int(p[P_NBALL]) patterns → p.geom.Z
    text = replace(text, "Int(p.geom.Z)" => "p.geom.Z")
    text = replace(text, "Float64(p.geom.Z)" => "Float64(p.geom.Z)")

    # Fix pilot_is_inner comparisons
    text = replace(text, "(p.cage.pilot_is_inner ? 1.0 : 0.0) > 0.5" => "p.cage.pilot_is_inner")
    
    # Fix params_tuple unpacking in kernel
    text = replace(text, "p = params_tuple[1]" => "p = params_tuple[1]::ODEParams")
    
    if changed
        write(path, text)
        println("  ✓ Patched: $path")
    else
        println("  ○ No changes: $path")
    end
end

println("Patching consumer files...")
patch_file!("src/Dynamics/kernel.jl")
patch_file!("src/Physics/kinematics.jl")
patch_file!("src/Dynamics/field_output.jl")

# Fix driver.jl: SimResult.params type + params_tuple wrapping
text = read("src/Dynamics/driver.jl", String)
text = replace(text, "params::Vector{Float64}" => "params::ODEParams")
text = replace(text, "params_tuple = (params,)  # wrap for kernel" => "params_tuple = (params,)  # wrap ODEParams for kernel")
# Fix postprocess references
text = replace(text, "Int(result.params[P_NBALL])" => "result.params.geom.Z")
text = replace(text, "result.params[P_OMEGA_IR]" => "result.params.load.omega_ir")
text = replace(text, "result.params[P_D]" => "result.params.geom.D")
text = replace(text, "result.params[P_DM]" => "result.params.geom.d_m")
text = replace(text, "result.params[P_ALPHA0]" => "result.params.geom.alpha_0")
write("src/Dynamics/driver.jl", text)
println("  ✓ Patched: src/Dynamics/driver.jl")

# Fix postprocess.jl
text = read("src/IO/postprocess.jl", String)
text = replace(text, "Int(result.params[P_NBALL])" => "result.params.geom.Z")
text = replace(text, "result.params[P_OMEGA_IR]" => "result.params.load.omega_ir")
text = replace(text, "result.params[P_D]" => "result.params.geom.D")
text = replace(text, "result.params[P_DM]" => "result.params.geom.d_m")
text = replace(text, "result.params[P_ALPHA0]" => "result.params.geom.alpha_0")
write("src/IO/postprocess.jl", text)
println("  ✓ Patched: src/IO/postprocess.jl")

# Fix the KinematicContext struct p field type
text = read("src/Physics/kinematics.jl", String)
text = replace(text, "p::Vector{Float64}" => "p::ODEParams")
write("src/Physics/kinematics.jl", text)
println("  ✓ Patched KinematicContext p type")

println("\n✅ All consumer files patched!")
