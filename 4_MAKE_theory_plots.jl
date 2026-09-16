###################################################################
############ 4. THEORY PLOT (SEASONAL WAGES SIMULATION) ###########
# Writes figures/simulation_predictions.pdf
# Run from the repository root: julia --project=. 4_MAKE_theory_plots.jl

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DataFrames
using Statistics
using Distributions
using Random
using StatsBase
using Plots
using Plots.PlotMeasures
using CSV

include(joinpath(@__DIR__, "code", "functions", "utility.jl"))

# The theory ABM (abm_cleaned.jl) has its own submodule set; code/abm/submodules belongs to the sweep ABM (test_seasons.jl).
submodule_dir = joinpath(@__DIR__, "code", "abm", "theory_submodules")
for file in readdir(submodule_dir)
    include(joinpath(submodule_dir, file))
end

include(joinpath(@__DIR__, "code", "abm", "abm_cleaned.jl"))

mkpath(joinpath(@__DIR__, "figures"))
include(joinpath(@__DIR__, "code", "plotting", "base_seasonality_theory_plot.jl"))
