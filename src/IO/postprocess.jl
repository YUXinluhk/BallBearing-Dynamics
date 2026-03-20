# =====================================================================
# IO/postprocess.jl — Result extraction and dimensional conversion
#
# Port of postprocess.py core extraction routines.
# =====================================================================

"""
    extract_ball_data(result, scales, ball_idx) → NamedTuple

Extract dimensional data for one ball from SimResult.
"""
function extract_ball_data(result::SimResult, ball_idx::Int)
    s = result.scales
    Z = result.params.geom.Z
    n_t = length(result.t)

    x = zeros(n_t)
    r = zeros(n_t)
    θ = zeros(n_t)
    θ̇ = zeros(n_t)
    ωx = zeros(n_t); ωy = zeros(n_t); ωz = zeros(n_t)

    for i in 1:n_t
        u = result.u[i, :]
        bp = ball_pos_view(u, ball_idx, Z)
        bv = ball_vel_view(u, ball_idx, Z)

        x[i]  = dim_length(s, bp[1])
        r[i]  = dim_length(s, bp[2])
        θ[i]  = bp[3]           # angle [rad]
        θ̇[i] = dim_angvel(s, bv[3])
        ωx[i] = dim_angvel(s, bv[4])
        ωy[i] = dim_angvel(s, bv[5])
        ωz[i] = dim_angvel(s, bv[6])
    end

    return (t=result.t, x=x, r=r, θ=θ, θ̇=θ̇, ωx=ωx, ωy=ωy, ωz=ωz)
end

"""
    extract_cage_data(result) → NamedTuple
"""
function extract_cage_data(result::SimResult)
    s = result.scales
    Z = result.params.geom.Z
    n_t = length(result.t)

    x = zeros(n_t); y = zeros(n_t); z = zeros(n_t)
    θ = zeros(n_t); θ̇ = zeros(n_t)

    for i in 1:n_t
        u = result.u[i, :]
        cp = cage_pos_view(u, Z)
        cv = cage_vel_view(u, Z)
        x[i] = dim_length(s, cp[1])
        y[i] = dim_length(s, cp[2])
        z[i] = dim_length(s, cp[3])
        θ[i] = cp[4]
        θ̇[i] = dim_angvel(s, cv[4])
    end

    return (t=result.t, x=x, y=y, z=z, θ=θ, θ̇=θ̇)
end

"""
    extract_ir_data(result) → NamedTuple
"""
function extract_ir_data(result::SimResult)
    s = result.scales
    Z = result.params.geom.Z
    n_t = length(result.t)

    x = zeros(n_t); y = zeros(n_t); z = zeros(n_t)

    for i in 1:n_t
        u = result.u[i, :]
        ip = ir_pos_view(u, Z)
        x[i] = dim_length(s, ip[1])
        y[i] = dim_length(s, ip[2])
        z[i] = dim_length(s, ip[3])
    end

    return (t=result.t, x=x, y=y, z=z)
end

"""
    compute_field_outputs(result) → Matrix

Compute field output for all saved timesteps.
Returns matrix of size (n_time, Z × N_FIELD_PER_BALL).
"""
function compute_field_outputs(result::SimResult)
    Z = result.params.geom.Z
    n_t = length(result.t)
    n_fo = Z * N_FIELD_PER_BALL

    fo_matrix = zeros(n_t, n_fo)

    for i in 1:n_t
        t_star = nondim_time(result.scales, result.t[i])
        # Reconstruct ComponentArray from flat vector using stored axes
        u_ca = ComponentArray(result.u[i, :], result.ca_axes...)
        outputs = field_output_kernel(t_star, u_ca, result.params)
        fo_matrix[i, :] .= flatten_field_outputs(outputs)
    end

    return fo_matrix
end

"""
    cage_speed_ratio(result) → Vector{Float64}

Compute cage speed / kinematic theory speed vs time.
"""
function cage_speed_ratio(result::SimResult)
    cage = extract_cage_data(result)
    ω_ir = result.params.load.omega_ir * result.scales.W  # dimensional
    γ = result.params.geom.D / result.params.geom.d_m
    cos_α₀ = cos(result.params.geom.alpha_0)
    ω_cage_kin = 0.5 * ω_ir * (1 - γ * cos_α₀)

    return cage.θ̇ ./ ω_cage_kin
end
