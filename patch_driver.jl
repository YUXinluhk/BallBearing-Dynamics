text = read("src/Dynamics/driver.jl", String)

old_cb = """
function quaternion_renormalize_callback(Z::Int)
    function affect!(integrator)
        u = integrator.u
        @inbounds for j in 1:Z
            off = ball_pos_offset(j) + 3
            s = u[off]
            v1 = u[off+1]
            v2 = u[off+2]
            v3 = u[off+3]
            n = sqrt(s^2 + v1^2 + v2^2 + v3^2)
            if n > 1e-15
                u[off] = s / n
                u[off+1] = v1 / n
                u[off+2] = v2 / n
                u[off+3] = v3 / n
            end
        end
    end
    DiscreteCallback((u, t, integrator) -> true, affect!)
end
"""

new_cb = """
function quaternion_renormalize_callback(Z::Int)
    function affect!(integrator)
        u = integrator.u
        @inbounds for j in 1:Z
            bp = u.ball[j].pos
            s = bp.q0
            v1 = bp.q1
            v2 = bp.q2
            v3 = bp.q3
            n = sqrt(s^2 + v1^2 + v2^2 + v3^2)
            if n > 1e-15
                bp.q0 = s / n
                bp.q1 = v1 / n
                bp.q2 = v2 / n
                bp.q3 = v3 / n
            end
        end
    end
    DiscreteCallback((u, t, integrator) -> true, affect!)
end
"""
text = replace(text, old_cb => new_cb)
write("src/Dynamics/driver.jl", text)
