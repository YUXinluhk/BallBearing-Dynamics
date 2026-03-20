# =====================================================================
# Physics/cage_contact.jl — Ball/cage pocket + cage/race pilot forces
#
# ⚠️  DEPRECATED — These functions are NOT called by the main ODE kernel.
#     kernel.jl implements its own inline pocket/pilot contact model with:
#       • Linear penalty stiffness (K·δ) instead of Hertz (K·δ^1.5)
#       • smooth_hertz_delta + tanh regularization (ODE-solver friendly)
#       • 3D slip friction with Newton's 3rd law back-reaction
#     This file is retained as a reference implementation only.
#
# Original port of dynamics_numba.py::_ball_cage_force_nb, _cage_race_force_nb
# =====================================================================

"""
    ball_cage_force(ball_θ, ball_x, ball_r, ball_ẋ,
                    cage_x, cage_y, cage_z, cage_θ, cage_ẋ, cage_ẏ, cage_ż, cage_θ̇,
                    pocket_r, pocket_clr, n_pockets, c_damp, μ_pocket, ball_idx)
    → (F_ball_x, F_ball_tang, F_cage_x, F_cage_y, F_cage_z, M_cage_θ)

Hertzian ball/cage pocket contact: F = K_h · δ^{1.5}.
"""
function ball_cage_force(ball_θ::Float64, ball_x::Float64, ball_r::Float64,
                         ball_ẋ::Float64,
                         cage_x::Float64, cage_y::Float64, cage_z::Float64,
                         cage_θ::Float64,
                         cage_ẋ::Float64, cage_ẏ::Float64, cage_ż::Float64,
                         cage_θ̇::Float64,
                         pocket_r::Float64, pocket_clr::Float64,
                         n_pockets::Int, c_damp::Float64,
                         μ_pocket::Float64, ball_idx::Int)
    # Pocket center azimuthal position
    ψ_pocket = cage_θ + (ball_idx - 1) * 2π / n_pockets

    # Tangential gap: ball azimuth vs pocket center
    Δθ = ball_θ - ψ_pocket
    # Wrap to [-π, π]
    Δθ -= 2π * round(Δθ / (2π))

    gap_tang = ball_r * Δθ   # tangential displacement [m]
    penetration = abs(gap_tang) - pocket_clr

    # Axial gap
    sψ, cψ = sincos(ψ_pocket)
    cage_r_at_pocket = -cage_y * sψ + cage_z * cψ
    gap_axial = ball_x - cage_x

    # Tangential force (Hertzian)
    F_tang = 0.0
    F_axial = 0.0

    if penetration > 0.0
        # Hertz K_h = (4/3)·E*·√R_eff, with E*≈1e11, R_eff≈pocket_r
        K_h = 1e10 * sqrt(pocket_r)  # simplified Hertz stiffness
        F_hertz = K_h * penetration^1.5

        # Relative tangential velocity
        v_rel_tang = ball_r * (ball_θ > ψ_pocket ? 1.0 : -1.0) -
                     ball_r * cage_θ̇  # simplified
        F_damp = c_damp * v_rel_tang

        F_tang = -sign(gap_tang) * (F_hertz + abs(F_damp))

        # Friction in axial direction
        F_axial = -μ_pocket * abs(F_hertz) * sign(gap_axial + 1e-30)
    end

    # Project forces
    F_ball_tang = F_tang
    F_ball_x = F_axial

    # Reaction on cage (Newton's 3rd law)
    F_cage_y = -F_tang * (-sψ)
    F_cage_z = -F_tang * cψ
    F_cage_x = -F_axial
    M_cage_θ = -F_tang * ball_r  # torque about bearing axis

    return (F_ball_x, F_ball_tang, F_cage_x, F_cage_y, F_cage_z, M_cage_θ)
end

"""
    cage_race_force(cage_y, cage_z, cage_ẏ, cage_ż, cage_θ̇, race_ω,
                    pilot_clr, c_damp, μ_pilot, pilot_inner,
                    cage_ir, cage_or)
    → (F_y, F_z, drag_torque)

Hertzian pilot contact + short-bearing hydrodynamic model.
"""
function cage_race_force(cage_y::Float64, cage_z::Float64,
                         cage_ẏ::Float64, cage_ż::Float64,
                         cage_θ̇::Float64, race_ω::Float64,
                         pilot_clr::Float64, c_damp::Float64,
                         μ_pilot::Float64, pilot_inner::Bool,
                         cage_ir::Float64, cage_or::Float64)
    # Eccentricity
    ecc = sqrt(cage_y^2 + cage_z^2)

    if ecc < 1e-15
        return (0.0, 0.0, 0.0)
    end

    # Unit direction
    ey = cage_y / ecc
    ez = cage_z / ecc

    F_r = 0.0
    if ecc > pilot_clr
        δ = ecc - pilot_clr
        # Hertz cylinder-in-bore
        K_pilot = 1e9 * sqrt(pilot_inner ? cage_ir : cage_or)
        F_r = K_pilot * δ^1.5

        # Damping
        v_radial = cage_ẏ * ey + cage_ż * ez
        F_r += c_damp * v_radial
    end

    F_y = -F_r * ey
    F_z = -F_r * ez

    # Drag torque from pilot ring
    R_pilot = pilot_inner ? cage_ir : cage_or
    Δω = cage_θ̇ - (pilot_inner ? race_ω : 0.0)
    drag_torque = -μ_pilot * abs(F_r) * sign(Δω + 1e-30) * R_pilot

    return (F_y, F_z, drag_torque)
end
