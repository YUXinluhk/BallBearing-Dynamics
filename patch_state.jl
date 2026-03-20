using ComponentArrays
using Statistics: mean

# Generate patch for src/Dynamics/state.jl
text = read("src/Dynamics/state.jl", String)

# 1. Imports
text = replace(text, "using StaticArrays" => "using StaticArrays\nusing ComponentArrays")

# 2. View accessors (We keep the same names but swap implementation to ComponentArrays)
# Find the "-- Zero-copy views -- " section
old_views = \"\"\"
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

"View ball j internal α₂ states [α₂_i, α₂_o]"
@inline function ball_alpha2_view(u, j, Z)
    off = ball_alpha2_offset(j, Z)
    @view u[off:off+N_ALPHA2_PER_BALL-1]
end

"View thermal nodes [T_i, T_o, T_b, T_oil]"
@inline function thermal_view(u, Z)
    off = thermal_offset(Z)
    @view u[off:off+N_THERMAL-1]
end

# ── Quaternion write-back ────────────────────────────────────────────

"Set ball j quaternion in state vector"
@inline function set_ball_quat!(u, j, Z, q::Quaternion)
    off = ball_pos_offset(j) + 3
    u[off] = real(q)
    u[off+1], u[off+2], u[off+3] = imag_part(q)
    return nothing
end
\"\"\"

new_views = \"\"\"
# ComponentArrays zero-copy views
@inline ir_pos_view(u, Z) = u.ir.pos
@inline ball_pos_view(u, j, Z) = u.ball[j].pos
@inline ball_trans_view(u, j, Z) = @view u.ball[j].pos[1:3]
@inline function ball_quat(u, j, Z)
    bp = u.ball[j].pos
    Quaternion{eltype(u)}(bp.q0, bp.q1, bp.q2, bp.q3)
end
@inline cage_pos_view(u, Z) = u.cage.pos

@inline ir_vel_view(u, Z) = u.ir.vel
@inline ball_vel_view(u, j, Z) = u.ball[j].vel
@inline function ball_omega(u, j, Z)
    bv = u.ball[j].vel
    SVector{3,eltype(u)}(bv.ωx, bv.ωy, bv.ωz)
end
@inline cage_vel_view(u, Z) = u.cage.vel

@inline ball_alpha2_view(u, j, Z) = u.ball[j].alpha2
@inline thermal_view(u, Z) = u.thermal

# ── Quaternion write-back ────────────────────────────────────────────
@inline function set_ball_quat!(u, j, Z, q::Quaternion)
    u.ball[j].pos.q0 = real(q)
    i1, i2, i3 = imag_part(q)
    u.ball[j].pos.q1 = i1
    u.ball[j].pos.q2 = i2
    u.ball[j].pos.q3 = i3
    return nothing
end
\"\"\"
text = replace(text, old_views => new_views)


# 3. Initialization
old_init = \"\"\"
function init_state(geom::BearingGeometry, qs::QuasiStaticResult,
    config::SimulationConfig, lub::LubricantParams)
    Z = geom.n_balls
    u0 = zeros(n_state(Z))

    # Inner race: axial displacement from QS
    ir = ir_pos_view(u0, Z)
    ir[1] = qs.delta_a  # x
    ir[2] = qs.delta_ry # y
    ir[3] = qs.delta_rz # z

    # Balls — place using averaged inner/outer contact angles
    # (QS solver 'Stalled' convergence → geometric inconsistency → average to
    # ensure both contacts have positive penetration at t=0)
    D = geom.d
    f_i = geom.f_i
    f_o = geom.f_o
    α₀ = alpha_free(geom)  # Use clearance-corrected free angle
    B_i = (f_i - 0.5) * D
    B_o = (f_o - 0.5) * D

    # Groove curvature center positions (fixed in bearing frame)
    x_ci = B_i * sin(α₀)                            # inner GC axial
    r_ci = geom.d_m / 2 + B_i * cos(α₀)             # inner GC radial
    x_co = -B_o * sin(α₀)                            # outer GC axial
    r_co = geom.d_m / 2 - B_o * cos(α₀)              # outer GC radial

    for j in 1:Z
        bp = ball_pos_view(u0, j, Z)
        ψ = (j - 1) * ball_spacing(geom)
        sψ, cψ = sincos(ψ)

        # Position from outer contact geometry
        α_o_j = deg2rad(qs.alpha_outer[j])
        L_o_j = B_o + qs.delta_outer[j]
        x_from_o = x_co + L_o_j * sin(α_o_j)
        r_from_o = r_co + L_o_j * cos(α_o_j)

        # Position from inner contact geometry (inner GC moves with IR)
        α_i_j = deg2rad(qs.alpha_inner[j])
        L_i_j = B_i + qs.delta_inner[j]
        r_ci_local = r_ci - qs.delta_ry * sψ + qs.delta_rz * cψ
        x_from_i = (x_ci + qs.delta_a) - L_i_j * sin(α_i_j)
        r_from_i = r_ci_local - L_i_j * cos(α_i_j)

        # Average to split geometric error between inner and outer contacts
        bp[1] = 0.5 * (x_from_o + x_from_i)          # x (axial)
        bp[2] = 0.5 * (r_from_o + r_from_i)          # r (radial)
        bp[3] = (j - 1) * ball_spacing(geom)          # θ (azimuth)

        # Quaternion: identity initially
        set_ball_quat!(u0, j, Z, QUAT_IDENTITY)

        # Ball velocity: orbital speed from QS
        bv = ball_vel_view(u0, j, Z)
        bv[3] = qs.omega_m[j]  # θ̇ = orbital speed
        # 【核心修复1】：用拟静力学节距角投影自旋角速度，彻底消除绝对侧翻导致的巨大初态冲击！
        beta_j = deg2rad(qs.beta[j])
        bv[4] = qs.omega_R[j] * cos(beta_j)  # ω_x (Roll轴向分量) 【Bug4修复】
        bv[5] = qs.omega_R[j] * sin(beta_j)  # ω_r (Pitch径向分量) 【Bug4修复】
        bv[6] = 0.0                          # ω_θ 绝对不可有宏观侧翻！
    end

    # Cage: at pitch radius, speed ≈ mean orbital
    cv = cage_vel_view(u0, Z)
    cv[4] = mean(qs.omega_m)  # θ̇_cage ≈ mean ball orbit

    # Thermal Nodes: array initialized to zero above, so we set to T_init
    tv = thermal_view(u0, Z)
    tv .= config.thermal.T_init

    # Heat accumulators: start at zero (∫H dt from t=0)
    # u0 is already zeros, so no explicit init needed.

    return u0
end
\"\"\"

new_init = \"\"\"
function init_state(geom::BearingGeometry, qs::QuasiStaticResult,
    config::SimulationConfig, lub::LubricantParams)
    Z = geom.n_balls

    # ComponentArrays requires vectors/tuples that are homogenous, or we define it explicitely
    ball_states = map(1:Z) do j
        (pos=(x=0.0, r=0.0, θ=0.0, q0=1.0, q1=0.0, q2=0.0, q3=0.0),
         vel=(x=0.0, r=0.0, θ=0.0, ωx=0.0, ωy=0.0, ωz=0.0),
         alpha2=(i=0.0, o=0.0))
    end
    
    u0_namedtuple = (
        ir = (pos=(x=qs.delta_a, y=qs.delta_ry, z=qs.delta_rz, γy=0.0, γz=0.0), vel=(x=0.0, y=0.0, z=0.0, γy=0.0, γz=0.0)),
        cage = (pos=(x=0.0, y=0.0, z=0.0, θ=0.0), vel=(x=0.0, y=0.0, z=0.0, θ=mean(qs.omega_m))),
        ball = ball_states,
        thermal = (Ti=config.thermal.T_init, To=config.thermal.T_init, Tb=config.thermal.T_init, Toil=config.thermal.T_init),
        heat = (ir=0.0, or_=0.0, b=0.0, oil=0.0, Q_amb=0.0, Cp_dT=0.0)
    )
    
    u0 = ComponentArray(u0_namedtuple)

    D = geom.d
    f_i = geom.f_i
    f_o = geom.f_o
    α₀ = alpha_free(geom)  # Use clearance-corrected free angle
    B_i = (f_i - 0.5) * D
    B_o = (f_o - 0.5) * D

    # Groove curvature center positions (fixed in bearing frame)
    x_ci = B_i * sin(α₀)                            # inner GC axial
    r_ci = geom.d_m / 2 + B_i * cos(α₀)             # inner GC radial
    x_co = -B_o * sin(α₀)                            # outer GC axial
    r_co = geom.d_m / 2 - B_o * cos(α₀)              # outer GC radial

    for j in 1:Z
        ψ = (j - 1) * ball_spacing(geom)
        sψ, cψ = sincos(ψ)

        # Position from outer contact geometry
        α_o_j = deg2rad(qs.alpha_outer[j])
        L_o_j = B_o + qs.delta_outer[j]
        x_from_o = x_co + L_o_j * sin(α_o_j)
        r_from_o = r_co + L_o_j * cos(α_o_j)

        # Position from inner contact geometry (inner GC moves with IR)
        α_i_j = deg2rad(qs.alpha_inner[j])
        L_i_j = B_i + qs.delta_inner[j]
        r_ci_local = r_ci - qs.delta_ry * sψ + qs.delta_rz * cψ
        x_from_i = (x_ci + qs.delta_a) - L_i_j * sin(α_i_j)
        r_from_i = r_ci_local - L_i_j * cos(α_i_j)

        bp = u0.ball[j].pos
        bp.x = 0.5 * (x_from_o + x_from_i)          # x (axial)
        bp.r = 0.5 * (r_from_o + r_from_i)          # r (radial)
        bp.θ = (j - 1) * ball_spacing(geom)          # θ (azimuth)
        
        # Quaternion is already identity in construction
        
        bv = u0.ball[j].vel
        bv.θ = qs.omega_m[j]  # θ̇ = orbital speed
        beta_j = deg2rad(qs.beta[j])
        bv.ωx = qs.omega_R[j] * cos(beta_j)  # ω_x (Roll轴向分量)
        bv.ωy = qs.omega_R[j] * sin(beta_j)  # ω_r (Pitch径向分量)
        bv.ωz = 0.0                          # ω_θ
    end

    return u0
end
\"\"\"
text = replace(text, old_init => new_init)

open("src/Dynamics/state.jl", "w") do io
    write(io, text)
end
println("Pacthed state.jl.")
