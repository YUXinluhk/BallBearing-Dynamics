# =====================================================================
# Dynamics/state.jl — SoA state vector layout with quaternion DOFs
#
# Position: IR(5) + Ball(7)×Z + Cage(4) = 9 + 7Z
# Velocity: IR(5) + Ball(6)×Z + Cage(4) = 9 + 6Z
# Total state = 18 + 13Z  (for Z=16 → 226)
#
# Ball position: [x, r, θ, q₀, q₁, q₂, q₃]  — 7 components
# Ball velocity: [ẋ, ṙ, θ̇, ω_x, ω_y, ω_z]  — 6 components
# =====================================================================

# ── Layout constants ─────────────────────────────────────────────────

const N_IR_POS = 5     # inner race: x, y, z, γ_y, γ_z
const N_BALL_POS = 7   # ball: x, r, θ, q₀, q₁, q₂, q₃
const N_CAGE_POS = 4   # cage: x, y, z, θ_cage
const N_IR_VEL = 5     # inner race: ẋ, ẏ, ż, γ̇_y, γ̇_z
const N_BALL_VEL = 6   # ball: ẋ, ṙ, θ̇, ω_x, ω_y, ω_z
const N_CAGE_VEL = 4   # cage: ẋ, ẏ, ż, θ̇_cage

"Total number of position DOFs"
n_pos_dofs(Z::Int) = N_IR_POS + N_BALL_POS * Z + N_CAGE_POS

"Total number of velocity DOFs"
n_vel_dofs(Z::Int) = N_IR_VEL + N_BALL_VEL * Z + N_CAGE_VEL

"Total state vector length"
n_state(Z::Int) = n_pos_dofs(Z) + n_vel_dofs(Z)

# ── Offset calculators ───────────────────────────────────────────────

# Position offsets (1-indexed)
@inline ir_pos_offset() = 1
@inline ball_pos_offset(j::Int) = N_IR_POS + (j - 1) * N_BALL_POS + 1
@inline cage_pos_offset(Z::Int) = N_IR_POS + N_BALL_POS * Z + 1

# Velocity offsets (1-indexed, after all positions)
@inline vel_base(Z::Int) = n_pos_dofs(Z)
@inline ir_vel_offset(Z::Int) = vel_base(Z) + 1
@inline ball_vel_offset(j::Int, Z::Int) = vel_base(Z) + N_IR_VEL + (j - 1) * N_BALL_VEL + 1
@inline cage_vel_offset(Z::Int) = vel_base(Z) + N_IR_VEL + N_BALL_VEL * Z + 1

# ── Zero-copy views ──────────────────────────────────────────────────

"View inner race position [x, y, z, γ_y, γ_z]"
@inline ir_pos_view(u, Z) = @view u[ir_pos_offset():ir_pos_offset()+N_IR_POS-1]

"View ball j position [x, r, θ, q₀, q₁, q₂, q₃]"
@inline function ball_pos_view(u, j, Z)
    off = ball_pos_offset(j)
    @view u[off:off+N_BALL_POS-1]
end

"View ball j translational position [x, r, θ]"
@inline function ball_trans_view(u, j, Z)
    off = ball_pos_offset(j)
    @view u[off:off+2]
end

"Extract ball j quaternion as Quaternion"
@inline function ball_quat(u, j, Z)
    off = ball_pos_offset(j) + 3  # skip x, r, θ
    Quaternion{eltype(u)}(u[off], u[off+1], u[off+2], u[off+3])
end

"View cage position [x, y, z, θ_cage]"
@inline function cage_pos_view(u, Z)
    off = cage_pos_offset(Z)
    @view u[off:off+N_CAGE_POS-1]
end

"View inner race velocity [ẋ, ẏ, ż, γ̇_y, γ̇_z]"
@inline ir_vel_view(u, Z) = @view u[ir_vel_offset(Z):ir_vel_offset(Z)+N_IR_VEL-1]

"View ball j velocity [ẋ, ṙ, θ̇, ω_x, ω_y, ω_z]"
@inline function ball_vel_view(u, j, Z)
    off = ball_vel_offset(j, Z)
    @view u[off:off+N_BALL_VEL-1]
end

"Ball j angular velocity [ω_x, ω_y, ω_z] as SVector"
@inline function ball_omega(u, j, Z)
    off = ball_vel_offset(j, Z) + 3  # skip ẋ, ṙ, θ̇
    SVector{3,eltype(u)}(u[off], u[off+1], u[off+2])
end

"View cage velocity [ẋ, ẏ, ż, θ̇_cage]"
@inline function cage_vel_view(u, Z)
    off = cage_vel_offset(Z)
    @view u[off:off+N_CAGE_VEL-1]
end

# ── Quaternion write-back ────────────────────────────────────────────

"Set ball j quaternion in state vector"
@inline function set_ball_quat!(u, j, Z, q::QuaternionF64)
    off = ball_pos_offset(j) + 3
    u[off] = real(q)
    u[off+1], u[off+2], u[off+3] = imag_part(q)
    return nothing
end

# ── State initialization ────────────────────────────────────────────

"""
    init_state(geom, qs_result, config) → u0

Build initial state vector from quasi-static solution.
"""
function init_state(geom::BearingGeometry, qs::QuasiStaticResult,
    config::SimulationConfig)
    Z = geom.n_balls
    u0 = zeros(n_state(Z))

    # Inner race: axial displacement from QS
    ir = ir_pos_view(u0, Z)
    ir[1] = qs.delta_a  # x

    # Balls — place using averaged inner/outer contact angles
    # (QS solver 'Stalled' convergence → geometric inconsistency → average to
    # ensure both contacts have positive penetration at t=0)
    D = geom.d
    f_i = geom.f_i
    f_o = geom.f_o
    α₀ = geom.alpha_0
    B_i = (f_i - 0.5) * D
    B_o = (f_o - 0.5) * D

    # Groove curvature center positions (fixed in bearing frame)
    x_ci = B_i * sin(α₀)                            # inner GC axial
    r_ci = geom.d_m / 2 + B_i * cos(α₀)             # inner GC radial
    x_co = -B_o * sin(α₀)                            # outer GC axial
    r_co = geom.d_m / 2 - B_o * cos(α₀)              # outer GC radial

    for j in 1:Z
        bp = ball_pos_view(u0, j, Z)

        # Position from outer contact geometry
        α_o_j = deg2rad(qs.alpha_outer[j])
        L_o_j = B_o + qs.delta_outer[j]
        x_from_o = x_co + L_o_j * sin(α_o_j)
        r_from_o = r_co + L_o_j * cos(α_o_j)

        # Position from inner contact geometry (inner GC moves with IR)
        α_i_j = deg2rad(qs.alpha_inner[j])
        L_i_j = B_i + qs.delta_inner[j]
        x_from_i = (x_ci + qs.delta_a) - L_i_j * sin(α_i_j)
        r_from_i = r_ci - L_i_j * cos(α_i_j)

        # Average to split geometric error between inner and outer contacts
        bp[1] = 0.5 * (x_from_o + x_from_i)          # x (axial)
        bp[2] = 0.5 * (r_from_o + r_from_i)          # r (radial)
        bp[3] = (j - 1) * ball_spacing(geom)          # θ (azimuth)

        # Quaternion: identity initially
        set_ball_quat!(u0, j, Z, QUAT_IDENTITY)

        # Ball velocity: orbital speed from QS
        bv = ball_vel_view(u0, j, Z)
        bv[3] = qs.omega_m[j]  # θ̇ = orbital speed
        # Spin: set ω_z ≈ spin speed about ball spin axis
        bv[6] = qs.omega_R[j]
    end

    # Cage: at pitch radius, speed ≈ mean orbital
    cv = cage_vel_view(u0, Z)
    cv[4] = mean(qs.omega_m)  # θ̇_cage ≈ mean ball orbit

    return u0
end
