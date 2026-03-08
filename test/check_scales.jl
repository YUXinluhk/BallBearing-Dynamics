using ADORE
geom, mat = bearing_7210B()
h_i, h_o = ADORE.create_bearing_hertz(geom, mat)
s = ADORE.Scales(geom, h_i)
cage = cage_from_bearing(geom)

println("Julia Scales (after T fix):")
println("  Q = ", s.Q)
println("  L = ", s.L)
println("  T = ", s.T)
println("  V = ", s.V)
println("  W = ", s.W)
println("  K = ", s.K)

m_ball = ADORE.ball_mass(geom)
m_ir = ADORE.inner_race_mass(geom)
m_cage = cage.cage_mass

println("\nNondim masses:")
println("  m_ball★ = ", ADORE.nondim_mass(s, m_ball))
println("  m_ir★   = ", ADORE.nondim_mass(s, m_ir))
println("  m_cage★ = ", ADORE.nondim_mass(s, m_cage))
println("  Y_i★    = ", ADORE.nondim_stiffness(h_i, s))
println("  Y_o★    = ", ADORE.nondim_stiffness(h_o, s))

# Compare with Python values
println("\nPython reference:")
println("  m_ball★ = 1.0")
println("  m_ir★   = 21.219")
println("  m_cage★ = 3.472")
println("  K_py    = Q/L = ", s.Q / s.L, " (Python stiffness scale)")
println("  K_jl    = Y   = ", s.K, " (Julia stiffness scale)")
