using ADORE
trac = traction_params_default()
println("Traction params: A=$(trac.A), B=$(trac.B), C=$(trac.C), D=$(trac.D)")
println()
for u in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 16.0, 50.0]
    k = ADORE.traction_coefficient(u, trac.A, trac.B, trac.C, trac.D)
    println("u=$(lpad(u, 6)) m/s  =>  kappa=$k")
end
