using ADORE

text = read("inputs/NASA_35mm_28k.toml", String)
text = replace(text, "t_end = 0.05" => "t_end = 1e-6")
text = replace(text, "dt_output = 1e-4" => "dt_output = 1e-6")
write("micro.toml", text)

geom, mat, lub, trac, cage, config = load_simulation_config("micro.toml")
sim = ADORE.run_simulation(geom, mat, lub, trac, cage, config)

println("SOLVED OK!")
