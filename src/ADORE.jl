"""
ADORE.jl — Angular-contact Deep-groove ball bearing Orbital Rolling Element dynamics.

Julia refactoring of ADORE_V2 Python code.
4-layer DOD architecture:
  1. Config/        — OOP types, validation
  2. Physics/       — Stateless pure functions
  3. QuasiStatic/   — 4Z+5 root-finding
  4. Dynamics/      — SoA state, ODE kernel
  + Transforms/     — Quaternion kinematics, frames
  + Scales/         — Nondimensionalization
  + IO/             — HDF5 persistence, postprocess
"""
module ADORE

using LinearAlgebra
using StaticArrays
using ComponentArrays
using Printf
using Statistics: mean
using SparseArrays

# ── Layer 1: Configuration ──────────────────────────────────────────
include("Config/validation.jl")
include("Config/types.jl")
include("Config/parser.jl")

export BearingGeometry, MaterialParams, LubricantParams, TractionParams
export CageGeometry, ChurningParams, IntegratorConfig, SimulationConfig
export ThermalParams, LPTNState
export r_ball, ball_mass, ball_inertia, ball_spacing, D_i, D_o
export inner_race_mass, sum_rho_inner, sum_rho_outer, F_rho_inner, F_rho_outer
export E_star, composite_modulus

# Convenience constructors
export bearing_7210B, bearing_7008C, bearing_7010C, bearing_7014C, bearing_custom,
  lubricant_mil_l_23699, traction_params_default
export cage_from_bearing, load_simulation_config

# ── Layer 2: Physics (stateless pure functions) ─────────────────────
# NOTE: ODEParams type definitions must be available before Physics/kinematics.jl
include("Dynamics/params_types.jl")
include("Physics/hertz.jl")
include("Physics/ehl.jl")
include("Physics/traction.jl")
include("Physics/drag_churning.jl")
include("Physics/cage_contact.jl")
include("Physics/rolling_resistance.jl")
include("Physics/thermal.jl")
include("Physics/kinematics.jl")

export hertz_F_of_k, solve_kappa_bisection
export hertz_Y_and_ab, hertz_from_curvatures, hertz_from_curvatures_vec
export hertz_contact_load, hertz_ab, smooth_hertz_delta
export HertzContact, create_bearing_hertz, nondim_stiffness
export viscosity, thermal_correction, film_thickness_hd
export traction_coefficient, spin_moment
export drag_coefficient_SN, drag_force, churning_moment, effective_density
export palmgren_Mrr, palmgren_power
export setup_thermal_capacitances, lptn_step, μ₀_at_temp, δr_thermal

# ── Layer 3: Transforms ─────────────────────────────────────────────
include("Transforms/quaternion.jl")
include("Transforms/frames.jl")

export QUAT_IDENTITY, quat_derivative, rotate_vector, inv_rotate_vector
export quat_to_rotmat, quat_renormalize, quat_from_euler_zyx
export omega_body_from_euler, quat_from_components, kinematics_quat_derivative_baumgarte
export T_azimuth, T_ac, contact_basis_from_displacements
export inertial_from_cylindrical, velocity_inertial_from_cylindrical

# ── Layer 4: Scales ─────────────────────────────────────────────────
include("Scales/nondim.jl")

export Scales
export nondim_force, nondim_length, nondim_time, nondim_vel, nondim_angvel
export nondim_moment, nondim_inertia, nondim_stiff, nondim_mass
export dim_force, dim_length, dim_time, dim_vel, dim_angvel, dim_moment, dim_inertia
export nondim_state

# ── Layer 5: Quasi-Static ───────────────────────────────────────────
include("QuasiStatic/result.jl")
include("QuasiStatic/solver.jl")

export QuasiStaticResult, solve_quasi_static

# ── Layer 6: Dynamics ───────────────────────────────────────────────
include("Dynamics/state.jl")
include("Dynamics/params.jl")
include("Dynamics/kernel.jl")
include("Dynamics/field_output.jl")
include("Dynamics/jacobian_sparsity.jl")
include("Dynamics/driver.jl")

export n_state, n_pos_dofs, n_vel_dofs
export ir_pos_view, ball_pos_view, ball_quat, ball_omega, cage_pos_view
export ir_vel_view, ball_vel_view, cage_vel_view, set_ball_quat!
export init_state, ode_rhs!
export build_params, ODEParams, N_PARAMS
export field_output_kernel, BallFieldOutput, N_FIELD_PER_BALL, flatten_field_outputs
export build_jacobian_sparsity
export SimResult, run_simulation, quaternion_renormalize_callback

# ── Layer 7: I/O ────────────────────────────────────────────────────
include("IO/buffer.jl")
include("IO/postprocess.jl")

export OutputBuffer, push_frame!, flush_buffer!, finalize_buffer!
export extract_ball_data, extract_cage_data, extract_ir_data
export compute_field_outputs, cage_speed_ratio

end # module ADORE
