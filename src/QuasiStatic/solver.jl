# =====================================================================
# QuasiStatic/solver.jl — 4Z+5 DOF Jones-Harris equilibrium solver
#
# Port of quasi_static.py.  Uses NonlinearSolve.jl (LevenbergMarquardt).
# =====================================================================

using NonlinearSolve

# ── Helper: smooth signed power ──────────────────────────────────────

"""
    eng_power(b, e; eps=1e-9)

C∞-smooth one-sided power: max(0, b)^e with softplus transition.
"""
@inline function eng_power(b, e; eps=1e-9)
    sb = 0.5 * (b + sqrt(b * b + eps * eps))
    return sb^e
end

# ── Curvature computation ────────────────────────────────────────────

function curv_inner(cos_alp_i::AbstractVector{T}, D, d_m, f_i) where T
    n = length(cos_alp_i)
    rho = Matrix{T}(undef, n, 4)
    for j in 1:n
        γ = D * cos_alp_i[j] / d_m
        rho[j, 1] = 2.0 / D
        rho[j, 2] = 2.0 / D
        rho[j, 3] = (2.0 / D) * (γ / (1.0 - γ))
        rho[j, 4] = -1.0 / (f_i * D)
    end
    return rho
end

function curv_outer(cos_alp_o::AbstractVector{T}, D, d_m, f_o) where T
    n = length(cos_alp_o)
    rho = Matrix{T}(undef, n, 4)
    for j in 1:n
        γ = D * cos_alp_o[j] / d_m
        rho[j, 1] = 2.0 / D
        rho[j, 2] = 2.0 / D
        rho[j, 3] = -(2.0 / D) * (γ / (1.0 + γ))
        rho[j, 4] = -1.0 / (f_o * D)
    end
    return rho
end

# ── Hertz stress (calc_stress) ───────────────────────────────────────

function calc_stress(delta::AbstractVector, rho::AbstractMatrix, E_star_val::Real)
    Z = length(delta)
    sum_rho = vec(sum(rho, dims=2))
    # Use Float64 curvatures for the bisection solver (not AD-dependent)
    sum_rho_f = Float64.(sum_rho)
    F_rho_arr_f = Float64.(vec(abs.((rho[:, 2] .+ rho[:, 4]) .- (rho[:, 1] .+ rho[:, 3])) ./ sum_rho))

    h = hertz_from_curvatures_vec(sum_rho_f, F_rho_arr_f, Float64(E_star_val))

    one_over_Et = 1.0 / E_star_val
    Q = h.Upsilon .* [eng_power(delta[j], 1.5) for j in 1:Z]

    base_args = 1.5 .* max.(Q, 0.0) ./ sum_rho .* one_over_Et
    base = [eng_power(base_args[j], 1.0 / 3.0) for j in 1:Z]
    a = h.a_star .* base
    b = h.b_star .* base
    ab = max.(a .* b, 1e-30)
    sigma = [ab[j] > 1e-30 ? 3.0 * Q[j] / (2π * ab[j]) : zero(Q[j]) for j in 1:Z]

    K_hertz = π .* h.kappa .* 2.0 ./ one_over_Et .*
              sqrt.(2.0 .* h.E2 ./ (9.0 .* sum_rho_f .* h.K1 .^ 3))

    return (Q=Q, sigma=sigma, a=a, E2=h.E2, K_hertz=K_hertz, b=b)
end

# ── Contact angle trigonometry ───────────────────────────────────────

function trig_fn(A_1, A_2, X_1, X_2, delta_i, delta_o, D, f_i, f_o)
    L_i = (f_i - 0.5) * D .+ delta_i
    L_o = (f_o - 0.5) * D .+ delta_o
    cos_alp_i = (A_2 .- X_2) ./ L_i
    cos_alp_o = X_2 ./ L_o
    sin_alp_i = (A_1 .- X_1) ./ L_i
    sin_alp_o = X_1 ./ L_o
    return (cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o)
end

# ── Kinematics ───────────────────────────────────────────────────────

function ball_orbit_speed(ω, sin_αi, cos_αi, sin_αo, cos_αo, γ, rc_id)
    cos_diff = cos_αi .* cos_αo .+ sin_αi .* sin_αo
    if rc_id == 0  # inner control
        return ω .* (cos_diff .- γ .* cos_αo) ./ (1.0 .+ cos_diff)
    else           # outer control
        return ω .* (1.0 .- γ .* cos_αi) ./ (1.0 .+ cos_diff)
    end
end

function gyro_moment(ω, sin_αi, cos_αi, sin_αo, cos_αo, γ, J_ball, rc_id)
    Z = length(sin_αi)
    if rc_id == 0
        β = atan.(sin_αi, cos_αi .- γ)
    else
        β = atan.(sin_αo, cos_αo .+ γ)
    end
    ω_m = ball_orbit_speed(ω, sin_αi, cos_αi, sin_αo, cos_αo, γ, rc_id)

    # Ball spin speed (simplified Jones-Harris race control)
    cos_β = cos.(β)
    tan_β = tan.(β)
    ω_R = similar(β)
    for j in 1:Z
        denom = ((cos_αo[j] + tan_β[j] * sin_αo[j]) / (1.0 + γ * cos_αo[j]) +
                 (cos_αi[j] + tan_β[j] * sin_αi[j]) / (1.0 - γ * cos_αi[j])) *
                γ * cos_β[j]
        denom = abs(denom) < 1e-30 ? 1e-30 : denom
        ω_R[j] = -ω / denom
    end

    M_g = J_ball .* ω_R .* ω_m .* sin.(β)
    return (M_g=M_g, ω_m=ω_m, ω_R=ω_R, β=β)
end

# ── Residual assembly ────────────────────────────────────────────────

"""
    QuasiStaticProblemParams — parameters passed to the residual function.
"""
struct QuasiStaticProblemParams
    Z::Int
    D::Float64
    d_m::Float64
    f_i::Float64
    f_o::Float64
    alpha_free::Float64
    gamma_tick::Float64
    R_i::Float64
    m_ball::Float64
    J_ball::Float64
    E_star::Float64
    psi_deg::Vector{Float64}
    # Load
    F_a::Float64
    F_ry::Float64
    F_rz::Float64
    Theta_y::Float64
    Theta_z::Float64
    omega::Float64
    delta_r_thermal::Float64
    rc_id::Int  # 0=inner, 1=outer
    lambda_i::Float64
    lambda_o::Float64
end

function quasi_static_residual!(F::AbstractVector, val_vec::AbstractVector,
    p::QuasiStaticProblemParams)
    Z = p.Z
    D, f_i, f_o = p.D, p.f_i, p.f_o
    cos_psi = cosd.(p.psi_deg)
    sin_psi = sind.(p.psi_deg)

    # ── Unpack state ──
    X_1 = val_vec[1:Z] .+ sin(p.alpha_free) * (f_o - 0.5) * D
    X_2 = val_vec[Z+1:2Z] .+ cos(p.alpha_free) * (f_o - 0.5) * D
    delta_i = val_vec[2Z+1:3Z]
    delta_o = val_vec[3Z+1:4Z]
    delta_a = val_vec[4Z+1]
    delta_ry = val_vec[4Z+2]
    delta_rz = val_vec[4Z+3]
    M_y = val_vec[4Z+4] * 1e4
    M_z = val_vec[4Z+5] * 1e4

    # ── Groove center distances ──
    R_i = p.R_i
    A_1 = (f_i + f_o - 1) * D * sin(p.alpha_free) .+ delta_a .+
          deg2rad(p.Theta_y) .* R_i .* cos_psi .+
          deg2rad(p.Theta_z) .* R_i .* sin_psi

    A_2 = (f_i + f_o - 1) * D * cos(p.alpha_free) .+
          delta_rz .* cos_psi .- delta_ry .* sin_psi .+ p.delta_r_thermal

    # ── Contact angles ──
    cos_αi, cos_αo, sin_αi, sin_αo = trig_fn(A_1, A_2, X_1, X_2,
        delta_i, delta_o, D, f_i, f_o)

    # ── Curvatures & Hertz ──
    ρ_ci = curv_inner(cos_αi, D, p.d_m, f_i)
    ρ_co = curv_outer(cos_αo, D, p.d_m, f_o)
    h_i = calc_stress(delta_i, ρ_ci, p.E_star)
    h_o = calc_stress(delta_o, ρ_co, p.E_star)

    # ── Kinematics ──
    gm = gyro_moment(p.omega, sin_αi, cos_αi, sin_αo, cos_αo,
        p.gamma_tick, p.J_ball, p.rc_id)
    M_g = gm.M_g
    ω_m = gm.ω_m
    F_c = 0.5 .* p.m_ball .* p.d_m .* ω_m .^ 2

    L_i = (f_i - 0.5) * D .+ delta_i
    L_o = (f_o - 0.5) * D .+ delta_o

    # Hertz loads via K_hertz
    Qi_15 = h_i.K_hertz .* [eng_power(delta_i[j], 1.5) for j in 1:Z]
    Qo_15 = h_o.K_hertz .* [eng_power(delta_o[j], 1.5) for j in 1:Z]

    λi, λo = p.lambda_i, p.lambda_o

    # ── Per-ball residuals ──
    # Scale geometric constraints by 1/(B²) to make dimensionless and comparable
    # to force residuals. Without scaling, Geo_err ~ O(1e-8 m²) appears small
    # but represents ~20% relative error → contact triangle doesn't close.
    B_i_sq = ((f_i - 0.5) * D)^2
    B_o_sq = ((f_o - 0.5) * D)^2
    geo_scale_i = 1.0 / B_i_sq   # ≈ 1.5e7 for 7210B
    geo_scale_o = 1.0 / B_o_sq
    Geo_i_err = ((A_1 .- X_1) .^ 2 .+ (A_2 .- X_2) .^ 2 .- L_i .^ 2) .* geo_scale_i
    Geo_o_err = (X_1 .^ 2 .+ X_2 .^ 2 .- L_o .^ 2) .* geo_scale_o

    F_b_ax = (λo .* M_g .* X_2 ./ D .- Qo_15 .* X_1) ./ L_o .+
             (Qi_15 .* (A_1 .- X_1) .- λi .* M_g ./ D .* (A_2 .- X_2)) ./ L_i

    F_b_rad = (Qo_15 .* X_2 .+ λo .* M_g .* X_1 ./ D) ./ L_o .-
              (Qi_15 .* (A_2 .- X_2) .+ λi .* M_g ./ D .* (A_1 .- X_1)) ./ L_i .- F_c

    # ── Global residuals (scaled to O(1) for LM residual balancing) ──
    Q_i_ax = (Qi_15 .* (A_1 .- X_1) .- λi .* M_g ./ D .* (A_2 .- X_2)) ./ L_i
    Q_i_rad = (Qi_15 .* (A_2 .- X_2) .+ λi .* M_g ./ D .* (A_1 .- X_1)) ./ L_i

    # Reference force scale: aligns force residuals ~O(1) with geometric residuals ~O(1)
    Q_ref = max(abs(p.F_a), maximum(abs, Qi_15), 1.0)
    M_ref = Q_ref * p.d_m  # moment reference scale

    F_ax_err = (p.F_a - sum(Q_i_ax)) / Q_ref
    F_rz_err = (p.F_rz - sum(Q_i_rad .* cos_psi)) / Q_ref
    F_ry_err = (p.F_ry + sum(Q_i_rad .* sin_psi)) / Q_ref
    M_by_err = (M_y - sum((Q_i_ax .* R_i .+ λi .* f_i .* M_g) .* cos_psi)) / M_ref
    M_bz_err = (M_z - sum((Q_i_ax .* R_i .+ λi .* f_i .* M_g) .* sin_psi)) / M_ref

    # ── Assemble (matching MATLAB ordering) ──
    F[1] = F_ax_err
    F[2] = F_ry_err
    F[3] = F_rz_err
    F[4] = M_by_err
    F[5] = M_bz_err
    F[6:5+Z] .= Geo_i_err
    F[5+Z+1:5+2Z] .= Geo_o_err
    F[5+2Z+1:5+3Z] .= F_b_ax ./ Q_ref
    F[5+3Z+1:5+4Z] .= F_b_rad ./ Q_ref

    return nothing
end

# ── Initial guess ────────────────────────────────────────────────────

function init_guess(Z::Int, args...)
    x0 = zeros(4Z + 5)
    x0[1:Z] .= 1e-5       # ΔX₁
    x0[Z+1:2Z] .= 1e-5    # ΔX₂
    x0[2Z+1:3Z] .= 1e-5   # δ_i
    x0[3Z+1:4Z] .= 1e-5   # δ_o
    x0[4Z+1] = 1e-5        # δ_a
    return x0
end

# ── Main solver ──────────────────────────────────────────────────────

"""
    solve_quasi_static(geom, mat; F_a, F_ry=0, F_rz=0, n_rpm=0, ...)
        → QuasiStaticResult

Solve (4Z+5) DOF Jones-Harris quasi-static equilibrium.
"""
function solve_quasi_static(geom::BearingGeometry, mat::MaterialParams;
    F_a::Float64=0.0, F_ry::Float64=0.0, F_rz::Float64=0.0,
    Theta_y::Float64=0.0, Theta_z::Float64=0.0,
    n_rpm::Float64=0.0,
    delta_r_thermal::Float64=0.0,
    x0::Union{Nothing,Vector{Float64}}=nothing,
    rc_assumption::Symbol=:outer,
    tol::Float64=1e-12,
    verbose::Bool=false)

    Z = geom.n_balls
    E_star = composite_modulus(mat, mat)
    γ = geom.d / geom.d_m
    R_i = 0.5 * geom.d_m + (geom.f_i - 0.5) * geom.d * cos(alpha_free(geom))
    ω = n_rpm / 60.0 * 2π
    psi = collect(0:Z-1) .* (360.0 / Z)

    rc_id = rc_assumption == :inner ? 0 : 1
    λi = rc_id == 0 ? 2.0 : 0.0
    λo = rc_id == 0 ? 0.0 : 2.0

    p = QuasiStaticProblemParams(
        Z, geom.d, geom.d_m, geom.f_i, geom.f_o, alpha_free(geom),
        γ, R_i, ball_mass(geom), ball_inertia(geom), E_star, psi,
        F_a, F_ry, F_rz, Theta_y, Theta_z, ω, delta_r_thermal,
        rc_id, λi, λo,
    )

    u0 = x0 === nothing ? init_guess(Z, alpha_free(geom), geom.d, geom.f_i, geom.f_o, F_a) : copy(x0)

    prob = NonlinearProblem(quasi_static_residual!, u0, p)
    sol = solve(prob, LevenbergMarquardt(; autodiff=AutoFiniteDiff());
        abstol=tol, reltol=tol, maxiters=5000)

    if verbose
        println("  NonlinearSolve success: $(sol.retcode)")
        println("  ||res||: $(norm(sol.resid))")
    end

    result = extract_qs_results(sol.u, p)
    # Accept Stalled as converged if residual is sufficiently small
    res_norm = norm(sol.resid)
    result.converged = (sol.retcode == ReturnCode.Success) ||
                       (sol.retcode == ReturnCode.Stalled && res_norm < 1e-4)
    result.sol_vec = copy(sol.u)

    return result
end

# ── Post-processing ──────────────────────────────────────────────────

function extract_qs_results(sol_vec::AbstractVector, p::QuasiStaticProblemParams)
    Z = p.Z
    D, f_i, f_o = p.D, p.f_i, p.f_o
    cos_psi = cosd.(p.psi_deg)
    sin_psi = sind.(p.psi_deg)

    X_1 = sol_vec[1:Z] .+ sin(p.alpha_free) * (f_o - 0.5) * D
    X_2 = sol_vec[Z+1:2Z] .+ cos(p.alpha_free) * (f_o - 0.5) * D
    delta_i = sol_vec[2Z+1:3Z]
    delta_o = sol_vec[3Z+1:4Z]
    delta_a = sol_vec[4Z+1]
    delta_ry = sol_vec[4Z+2]
    delta_rz = sol_vec[4Z+3]
    M_y_val = sol_vec[4Z+4] * 1e4
    M_z_val = sol_vec[4Z+5] * 1e4

    A_1 = (f_i + f_o - 1) * D * sin(p.alpha_free) .+ delta_a .+
          deg2rad(p.Theta_y) .* p.R_i .* cos_psi .+
          deg2rad(p.Theta_z) .* p.R_i .* sin_psi
    A_2 = (f_i + f_o - 1) * D * cos(p.alpha_free) .+
          delta_rz .* cos_psi .+ delta_ry .* sin_psi .+ p.delta_r_thermal

    cos_αi, cos_αo, sin_αi, sin_αo = trig_fn(A_1, A_2, X_1, X_2,
        delta_i, delta_o, D, f_i, f_o)

    ρ_ci = curv_inner(cos_αi, D, p.d_m, f_i)
    ρ_co = curv_outer(cos_αo, D, p.d_m, f_o)
    h_i = calc_stress(delta_i, ρ_ci, p.E_star)
    h_o = calc_stress(delta_o, ρ_co, p.E_star)

    gm = gyro_moment(p.omega, sin_αi, cos_αi, sin_αo, cos_αo,
        p.gamma_tick, p.J_ball, p.rc_id)

    F_c = 0.5 .* p.m_ball .* p.d_m .* gm.ω_m .^ 2

    result = QuasiStaticResult(Z)
    result.delta_inner = collect(delta_i)
    result.delta_outer = collect(delta_o)
    result.Q_inner = h_i.Q
    result.Q_outer = h_o.Q
    result.alpha_inner = rad2deg.(atan.(sin_αi, cos_αi))
    result.alpha_outer = rad2deg.(atan.(sin_αo, cos_αo))
    result.sigma_max_inner = h_i.sigma
    result.sigma_max_outer = h_o.sigma
    result.F_c = F_c
    result.M_g = gm.M_g
    result.delta_a = delta_a
    result.delta_ry = delta_ry
    result.delta_rz = delta_rz
    result.M_y = M_y_val
    result.M_z = M_z_val
    result.race_control = p.rc_id == 0 ? :inner : :outer
    result.omega_m = gm.ω_m
    result.omega_R = gm.ω_R
    result.beta = rad2deg.(gm.β)

    return result
end
