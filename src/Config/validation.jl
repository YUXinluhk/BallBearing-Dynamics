# =====================================================================
# Config/validation.jl — Boundary-check helpers
#
# Julia port of Python validation.py.
# Throws ArgumentError on invalid input; issues @warn for unusual ranges.
# =====================================================================

"""
    check_positive(name, val, ctx)

Throw `ArgumentError` if `val ≤ 0`.
"""
function check_positive(name::AbstractString, val::Real, ctx::AbstractString="")
    val > 0 && return nothing
    throw(ArgumentError("$name ($ctx) must be > 0, got $val"))
end

"""
    check_non_negative(name, val, ctx)

Throw `ArgumentError` if `val < 0`.
"""
function check_non_negative(name::AbstractString, val::Real, ctx::AbstractString="")
    val >= 0 && return nothing
    throw(ArgumentError("$name ($ctx) must be >= 0, got $val"))
end

"""
    check_range(name, val, lo, hi, ctx; inclusive=:both)

Throw `ArgumentError` if `val` is outside `[lo, hi]`.
`inclusive` ∈ {:both, :neither, :left, :right}.
"""
function check_range(name::AbstractString, val::Real, lo::Real, hi::Real,
                     ctx::AbstractString=""; inclusive::Symbol=:both)
    ok = if inclusive == :both
        lo <= val <= hi
    elseif inclusive == :neither
        lo < val < hi
    elseif inclusive == :left
        lo <= val < hi
    elseif inclusive == :right
        lo < val <= hi
    else
        error("Unknown inclusive mode: $inclusive")
    end
    ok && return nothing
    throw(ArgumentError("$name ($ctx) must be in ($lo, $hi), got $val"))
end

"""
    warn_range(name, val, lo, hi, ctx)

Issue `@warn` if `val` is outside the *typical* range `[lo, hi]`.
"""
function warn_range(name::AbstractString, val::Real, lo::Real, hi::Real,
                    ctx::AbstractString="")
    if val < lo || val > hi
        @warn "$name ($ctx) = $val is outside typical range [$lo, $hi]"
    end
    return nothing
end

"""
    validate_simulation(geom, mat; cage_geom=nothing, lub=nothing, trac=nothing, config=nothing)

Cross-object consistency checks.
"""
function validate_simulation(geom, mat; cage_geom=nothing, lub=nothing,
                             trac=nothing, config=nothing)
    # Cage pocket count must equal number of balls
    if cage_geom !== nothing && cage_geom.n_pockets != geom.n_balls
        throw(ArgumentError(
            "pocket count ($(cage_geom.n_pockets)) must equal n_balls ($(geom.n_balls))"))
    end
    return nothing
end
