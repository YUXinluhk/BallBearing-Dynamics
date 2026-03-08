# =====================================================================
# QuasiStatic/result.jl — Result type for quasi-static equilibrium
# =====================================================================

"""
    QuasiStaticResult — Full results from quasi-static equilibrium.

All per-ball arrays are length Z.
"""
mutable struct QuasiStaticResult
    converged::Bool
    n_iterations::Int
    exit_flag::Int

    # Per-ball arrays (length Z)
    delta_inner::Vector{Float64}
    delta_outer::Vector{Float64}
    Q_inner::Vector{Float64}
    Q_outer::Vector{Float64}
    alpha_inner::Vector{Float64}   # [deg]
    alpha_outer::Vector{Float64}   # [deg]
    sigma_max_inner::Vector{Float64}
    sigma_max_outer::Vector{Float64}
    F_c::Vector{Float64}
    M_g::Vector{Float64}

    # Global displacements
    delta_a::Float64
    delta_ry::Float64
    delta_rz::Float64
    M_y::Float64
    M_z::Float64

    # Kinematics
    race_control::Symbol   # :inner or :outer
    omega_m::Vector{Float64}
    omega_R::Vector{Float64}
    beta::Vector{Float64}          # [deg]
    sol_vec::Vector{Float64}       # raw solution [4Z+5] for continuation
end

"Create empty QuasiStaticResult for Z balls"
function QuasiStaticResult(Z::Int)
    QuasiStaticResult(
        false, 0, 0,
        zeros(Z), zeros(Z), zeros(Z), zeros(Z),
        zeros(Z), zeros(Z), zeros(Z), zeros(Z),
        zeros(Z), zeros(Z),
        0.0, 0.0, 0.0, 0.0, 0.0,
        :outer,
        zeros(Z), zeros(Z), zeros(Z),
        zeros(4Z + 5),
    )
end
