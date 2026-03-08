using Pkg
deps = Pkg.dependencies()
targets = ["Plots", "CairoMakie", "Makie", "PyPlot", "GR", "UnicodePlots", "PythonPlot"]
for (uuid, info) in deps
    if info.name in targets
        println(info.name)
    end
end
