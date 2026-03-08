# =====================================================================
# Scales/nondim.jl — Non-dimensionalization scales (Gupta Ch. 6)
#
# Port of nondim.py::Scales + all conversion functions.
# =====================================================================

"""
    Scales — characteristic scales for non-dimensionalization.

Q₀ = Y_inner·L₀^{3/2}  (force scale from inner Hertz stiffness)
L₀ = r_ball             (length scale)
T₀ = √(m_ball / Y_inner) · L₀^{-1/4}  (time scale)
"""
struct Scales
    Q::Float64    # force [N]
    L::Float64    # length [m]
    T::Float64    # time [s]
    V::Float64    # velocity [m/s]
    W::Float64    # angular velocity [rad/s]
    K::Float64    # stiffness [N/m^{3/2}]
    I::Float64    # moment of inertia [kg·m²]
    M::Float64    # moment [N·m]
end

"""
    Scales(geom, hertz_inner) → Scales

Create non-dimensional scales from bearing geometry + inner Hertz contact.
"""
function Scales(geom::BearingGeometry, hertz_inner::HertzContact)
    L = r_ball(geom)                    # length scale = ball radius
    Y = hertz_inner.Upsilon
    m = ball_mass(geom)

    Q = Y * L^1.5                       # force scale
    T = sqrt(m / (Y * sqrt(L)))          # time scale: √(m/(Y·L^{1/2})) = √(m/Y)·L^{-1/4}
    V = L / T                           # velocity scale
    W = 1.0 / T                         # angular velocity scale
    K = Y                               # stiffness scale
    I_scale = m * L^2                   # inertia scale
    M_scale = Q * L                     # moment scale

    Scales(Q, L, T, V, W, K, I_scale, M_scale)
end

# ── Non-dimensionalize ────────────────────────────────────────────────

@inline nondim_force(s::Scales, F) = F / s.Q
@inline nondim_length(s::Scales, x) = x / s.L
@inline nondim_time(s::Scales, t) = t / s.T
@inline nondim_vel(s::Scales, v) = v / s.V
@inline nondim_angvel(s::Scales, ω) = ω / s.W
@inline nondim_moment(s::Scales, M) = M / s.M
@inline nondim_inertia(s::Scales, I) = I / s.I
@inline nondim_stiff(s::Scales, K) = K / s.K
@inline nondim_mass(s::Scales, m) = m * s.L / (s.Q * s.T^2)

# ── Dimensionalize ────────────────────────────────────────────────────

@inline dim_force(s::Scales, F★) = F★ * s.Q
@inline dim_length(s::Scales, x★) = x★ * s.L
@inline dim_time(s::Scales, t★) = t★ * s.T
@inline dim_vel(s::Scales, v★) = v★ * s.V
@inline dim_angvel(s::Scales, ω★) = ω★ * s.W
@inline dim_moment(s::Scales, M★) = M★ * s.M
@inline dim_inertia(s::Scales, I★) = I★ * s.I

# ── State vector conversion ──────────────────────────────────────────

"""
    nondim_state(u_dim, scales, n_balls) → u_star

Convert full dimensional state vector to non-dimensional form.
"""
function nondim_state(u_dim::AbstractVector{Float64}, scales::Scales, n_balls::Int)
    u★ = similar(u_dim)
    # Layout: positions then velocities
    # Positions: lengths and angles
    # Inner race: 5 pos (x, y, z = lengths; γ_y, γ_z = angles [dimensionless])
    # Balls: 7 per ball (x, r, θ = length, length, angle; q₀,q₁,q₂,q₃ = unitless)
    # Cage: 4 pos (x, y, z = lengths; θ_cage = angle)

    N_ir = 5
    N_ball = 7
    N_cage = 4
    N_pos = N_ir + N_ball * n_balls + N_cage

    # Velocities:
    N_ir_v = 5
    N_ball_v = 6
    N_cage_v = 4

    s = scales

    # ── Positions ──
    # Inner race [x, y, z, γ_y, γ_z]
    for i in 1:3
        u★[i] = u_dim[i] / s.L          # lengths
    end
    u★[4] = u_dim[4]                     # angles = dimensionless
    u★[5] = u_dim[5]

    # Balls
    for j in 1:n_balls
        off = N_ir + (j - 1) * N_ball
        u★[off+1] = u_dim[off+1] / s.L  # x
        u★[off+2] = u_dim[off+2] / s.L  # r
        u★[off+3] = u_dim[off+3]         # θ (angle, dimensionless)
        for k in 4:7
            u★[off+k] = u_dim[off+k]     # quaternion components (dimensionless)
        end
    end

    # Cage [x, y, z, θ]
    off_c = N_ir + N_ball * n_balls
    u★[off_c+1] = u_dim[off_c+1] / s.L
    u★[off_c+2] = u_dim[off_c+2] / s.L
    u★[off_c+3] = u_dim[off_c+3] / s.L
    u★[off_c+4] = u_dim[off_c+4]         # angle

    # ── Velocities ──
    vel_start = N_pos  # 0-based offset for velocity section

    # IR velocities [ẋ, ẏ, ż (linear), γ̇_y, γ̇_z (angular)]
    for i in 1:3
        u★[vel_start+i] = u_dim[vel_start+i] / s.V
    end
    u★[vel_start+4] = u_dim[vel_start+4] / s.W
    u★[vel_start+5] = u_dim[vel_start+5] / s.W

    # Ball velocities [ẋ, ṙ, θ̇, ω_x, ω_y, ω_z]
    for j in 1:n_balls
        off_v = vel_start + N_ir_v + (j - 1) * N_ball_v
        u★[off_v+1] = u_dim[off_v+1] / s.V  # ẋ
        u★[off_v+2] = u_dim[off_v+2] / s.V  # ṙ
        u★[off_v+3] = u_dim[off_v+3] / s.W  # θ̇
        for k in 4:6
            u★[off_v+k] = u_dim[off_v+k] / s.W  # ω_x, ω_y, ω_z
        end
    end

    # Cage velocities [ẋ, ẏ, ż, θ̇]
    off_cv = vel_start + N_ir_v + N_ball_v * n_balls
    u★[off_cv+1] = u_dim[off_cv+1] / s.V
    u★[off_cv+2] = u_dim[off_cv+2] / s.V
    u★[off_cv+3] = u_dim[off_cv+3] / s.V
    u★[off_cv+4] = u_dim[off_cv+4] / s.W

    return u★
end
