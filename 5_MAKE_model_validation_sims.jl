###################################################################
############ 5. MODEL VALIDATION - SIMULATED SWEEPS ###############
# Writes data/sweeps/{earnings_on_kesi,weather_on_earnings}/abm/*.csv and
#        the two sweep_list.csv files read by 6_MAKE.
# KESI_WORKERS sets the number of worker processes (default 20; the full
# sweep takes about a day on 20 cores). KESI_SMOKE=1 runs a 2x2 grid of short
# simulations to check the pipeline.
# Run from the repository root: julia --project=. 5_MAKE_model_validation_sims.jl

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Distributed
const DEFAULT_WORKERS = 20
addprocs(parse(Int, get(ENV, "KESI_WORKERS", string(DEFAULT_WORKERS))); exeflags = "--project=$(@__DIR__)")

@everywhere const KESI_SMOKE = haskey(ENV, "KESI_SMOKE")
@everywhere const SMOKE_YEARS = 50
@everywhere smoke_grid(values) = KESI_SMOKE ? values[[1, end]] : values

@everywhere using DataFrames
@everywhere using Statistics
@everywhere using Distributions
@everywhere using Random
@everywhere using StatsBase
@everywhere using StatsFuns
@everywhere using CSV

@everywhere const REPO_ROOT = $(@__DIR__)
@everywhere cd(REPO_ROOT)

@everywhere include(joinpath(REPO_ROOT, "code", "functions", "utility.jl"))

@everywhere submodule_dir = joinpath(REPO_ROOT, "code", "abm", "submodules")
@everywhere for file in readdir(submodule_dir)
    include(joinpath(submodule_dir, file))
end

@everywhere include(joinpath(REPO_ROOT, "code", "abm", "test_seasons.jl"))

for sweep in ("earnings_on_kesi", "weather_on_earnings"), stage in ("abm", "stan")
    mkpath(joinpath(REPO_ROOT, "data", "sweeps", sweep, stage))
end

include(joinpath(REPO_ROOT, "code", "sweeps", "earnings_on_kesi_sweep.jl"))
include(joinpath(REPO_ROOT, "code", "sweeps", "weather_on_earnings_sweep.jl"))
