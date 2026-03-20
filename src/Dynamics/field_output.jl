# =====================================================================
# Dynamics/field_output.jl — Field output extraction (structured approach)
# 完全镜像同步 kernel 动力学内核 (deepthink.md 参考架构)
# =====================================================================

"""
    BallFieldOutput — Structured field output for a single ball.

All fields are dimensional (SI units) unless noted.
Replaces the legacy 42-constant flat array indexing.
"""
struct BallFieldOutput
    # Contact loads [N]
    Q_i::Float64;       Q_o::Float64
    # Contact angles [rad]
    α_i::Float64;       α_o::Float64
    # Contact deflections [m]
    δ_i::Float64;       δ_o::Float64
    # Contact semi-axes [m]
    a_i::Float64;       b_i::Float64
    a_o::Float64;       b_o::Float64
    # Sliding velocities [m/s]
    u_slide_i::Float64; u_slide_o::Float64
    # Mean entrainment velocities [m/s]
    u_mean_i::Float64;  u_mean_o::Float64
    # Film thicknesses [m]
    h_film_i::Float64;  h_film_o::Float64
    # Traction coefficients [-]
    κ_i::Float64;       κ_o::Float64
    # Traction forces [N]
    F_trac_i::Float64;  F_trac_o::Float64
    # Spin moments [N·m]
    M_spin_i::Float64;  M_spin_o::Float64
    # Drag force [N] and churning moment [N·m]
    F_drag::Float64;    M_churn::Float64
    # Ball spin components [rad/s]
    ω_x::Float64;       ω_y::Float64;       ω_z::Float64
    # Ball orbital position [m] and speed [rad/s]
    r_ball::Float64;    θ_dot::Float64
    # Ball quaternion
    q0::Float64; q1::Float64; q2::Float64; q3::Float64
    # Heat generation [W]
    H_slide_i::Float64; H_slide_o::Float64
    H_spin_i::Float64;  H_spin_o::Float64
    H_drag::Float64;    H_churn::Float64
    # Pocket force [N]
    F_pocket::Float64
    # Spin speeds [rad/s]
    ω_spin_i::Float64;  ω_spin_o::Float64
end

# Number of fields per ball (for backward compat with flat matrix)
const N_FIELD_PER_BALL = fieldcount(BallFieldOutput)

"""
    field_output_kernel(t, u, p) → Vector{BallFieldOutput}

Extract diagnostic field outputs from state vector u at time t.
完全镜像同步：与 kernel.jl 使用完全相同的接触几何和运动学公式。
"""
function field_output_kernel(t::Float64, u::AbstractVector, p::ODEParams)
    Z = p.geom.Z
    outputs = Vector{BallFieldOutput}(undef, Z)

    ctx = build_kinematic_context(u, p, t)

    @inbounds for j in 1:Z
        res = compute_ball_tehd_kinematics(ctx, j)

        bp = u.ball[j].pos
        bv = u.ball[j].vel
        x_b, r_b, θ_b = bp.x, bp.r, bp.θ
        ẋ_b, ṙ_b, θ̇_b = bv.x, bv.r, bv.θ
        q = ball_quat(u, j, Z)

        α_i = atan(res.dx_i, res.dr_i)
        α_o = atan(res.dx_o, res.dr_o)

        h_film_i = film_thickness_hd(res.u_entrain_i, res.Q_i_dim, ctx.E_prime_dim, res.R_ball_k * ctx.L_scale,
            (res.b_i_dim > 0 ? res.a_i_dim / res.b_i_dim : 1.0), ctx.μ_oil_dyn, p.lub.alpha_pv, p.lub.K_th, ctx.T_oil_node)
        h_film_o = film_thickness_hd(res.u_entrain_o, res.Q_o_dim, ctx.E_prime_dim, res.R_ball_k * ctx.L_scale,
            (res.b_o_dim > 0 ? res.a_o_dim / res.b_o_dim : 1.0), ctx.μ_oil_dyn, p.lub.alpha_pv, p.lub.K_th, ctx.T_oil_node)

        κ_i = res.F_trac_i_dim / (res.Q_i_dim + 1e-24)
        κ_o = res.F_trac_o_dim / (res.Q_o_dim + 1e-24)

        pocket_θ = ctx.θ_cg + (j - 1) * 2π / Z
        Δθ = θ_b - pocket_θ
        Δθ -= 2π * round(Δθ / (2π))
        gap = r_b * Δθ
        pen = hypot(gap, 1e-12) - p.cage.pocket_clr
        F_pk = pen > 0 ? p.cage.k_pocket * pen * ctx.Q_scale : 0.0

        H_drag = (res.P_drag_nd + res.P_ehl_roll_nd) * ctx.Q_scale * ctx.V_scale
        H_churn = res.P_drag_spin_nd * ctx.Q_scale * ctx.V_scale

        q1_v, q2_v, q3_v = imag_part(q)

        outputs[j] = BallFieldOutput(
            res.Q_i_dim, res.Q_o_dim,                   # Q_i, Q_o [N]
            α_i, α_o,                                     # α_i, α_o [rad]
            res.δ_i_sm * ctx.L_scale, res.δ_o_sm * ctx.L_scale,  # δ [m]
            res.a_i_dim, res.b_i_dim, res.a_o_dim, res.b_o_dim,  # semi-axes [m]
            res.u_roll_i_dim, res.u_roll_o_dim,           # u_slide [m/s]
            0.5 * (res.v_ir_θ_dim + res.v_ball_θ * ctx.V_scale),  # u_mean_i [m/s]
            0.5 * (res.v_or_θ_dim + res.v_ball_θ * ctx.V_scale),  # u_mean_o [m/s]
            h_film_i, h_film_o,                           # h_film [m]
            κ_i, κ_o,                                     # traction coeff [-]
            res.F_trac_i_dim, res.F_trac_o_dim,           # F_trac [N]
            res.M_sp_i_dim, res.M_sp_o_dim,               # M_spin [N·m]
            res.F_drag_dim, res.M_ch_dim,                 # drag/churn
            res.ω_ball_x_dim, res.ω_ball_y_dim, res.ω_ball_z_dim,  # ω [rad/s]
            r_b * ctx.L_scale, θ̇_b * ctx.W_scale,        # r, θ̇
            real(q), q1_v, q2_v, q3_v,                    # quaternion
            res.H_slide_i_dim, res.H_slide_o_dim,         # H_slide [W]
            res.H_spin_i_dim, res.H_spin_o_dim,           # H_spin [W]
            H_drag, H_churn,                              # H_drag/churn [W]
            F_pk,                                          # F_pocket [N]
            res.ω_spin_i_dim, res.ω_spin_o_dim            # ω_spin [rad/s]
        )
    end

    return outputs
end

# ── Flat vector conversion for backward compat with compute_field_outputs ──
"""
    flatten_field_outputs(outputs::Vector{BallFieldOutput}) → Vector{Float64}

Convert structured outputs to flat vector for backward-compatible matrix storage.
"""
function flatten_field_outputs(outputs::Vector{BallFieldOutput})
    Z = length(outputs)
    fo = zeros(Z * N_FIELD_PER_BALL)
    for j in 1:Z
        base = (j - 1) * N_FIELD_PER_BALL
        o = outputs[j]
        for (k, fname) in enumerate(fieldnames(BallFieldOutput))
            fo[base + k] = Float64(getfield(o, fname))
        end
    end
    return fo
end
