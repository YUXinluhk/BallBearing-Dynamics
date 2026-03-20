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

Convert full dimensional ComponentArray state vector to non-dimensional form.
Uses semantic property access — completely immune to state layout changes.
"""
function nondim_state(u_dim::AbstractVector{Float64}, scales::Scales, n_balls::Int)
    u★ = copy(u_dim)
    s = scales

    # ── Inner race position (lengths → nondim, angles pass through) ──
    u★.ir.pos.x /= s.L
    u★.ir.pos.y /= s.L
    u★.ir.pos.z /= s.L
    # γy, γz are angles — already dimensionless

    # ── Inner race velocity (linear → V, angular → W) ──
    u★.ir.vel.x /= s.V
    u★.ir.vel.y /= s.V
    u★.ir.vel.z /= s.V
    u★.ir.vel.γy /= s.W
    u★.ir.vel.γz /= s.W

    # ── Ball states ──
    for j in 1:n_balls
        # Position: x, r are lengths; θ and quaternion are dimensionless
        u★.ball[j].pos.x /= s.L
        u★.ball[j].pos.r /= s.L
        # θ, q0-q3 are dimensionless — pass through

        # Velocity: ẋ, ṙ are linear; θ̇, ωx, ωy, ωz are angular
        u★.ball[j].vel.x /= s.V
        u★.ball[j].vel.r /= s.V
        u★.ball[j].vel.θ /= s.W
        u★.ball[j].vel.ωx /= s.W
        u★.ball[j].vel.ωy /= s.W
        u★.ball[j].vel.ωz /= s.W
    end

    # ── Cage position (lengths → nondim, angle pass through) ──
    u★.cage.pos.x /= s.L
    u★.cage.pos.y /= s.L
    u★.cage.pos.z /= s.L
    # θ is angle — dimensionless

    # ── Cage velocity ──
    u★.cage.vel.x /= s.V
    u★.cage.vel.y /= s.V
    u★.cage.vel.z /= s.V
    u★.cage.vel.θ /= s.W

    # ── Thermal & heat accumulator: pass through unchanged ──
    # Temperatures [K] are dimensional, heat accumulators are nondim energy

    return u★
end
