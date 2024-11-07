####################
#LOAD PACKAGES #####

###############################################
########## USER SET NUMBER OF CORES ###########
addprocs(20)


######################################
############ DETERMINE COMPUTER ######
using Distributed

# Set up nprocs()
@everywhere using DataFrames
@everywhere using Statistics
@everywhere using Distributions
@everywhere using Random
@everywhere using Distributions
@everywhere using StatsBase
@everywhere using Plots
@everywhere using Plots.PlotMeasures
@everywhere using JLD2
@everywhere using Serialization
@everywhere using Statistics
@everywhere using ColorSchemes
@everywhere using GLM
@everywhere using StatsFuns
@everywhere using CSV
# SET working directory as per readme   

#####################################
######## Initalize Functions ########

@everywhere include(string(pwd(), "\\code\\functions\\utility.jl"))

######################################
#### Initalize submodules ############

@everywhere files = readdir(string(pwd(), ("\\code\\abm\\submodules")))
@everywhere for i in files  include(string(pwd(), "\\code\\abm\\submodules\\$i")) end

######################################
######### CHOOSE ABM VERSION #########
@everywhere include(string(pwd(), "\\code\\abm\\test_seasons.jl"))

########################################
######## Simulate data #################
@everywhere include("\\code\\sweeps\\earnings_on_kesi.jl")
@everywhere include("\\code\\sweeps\\weather_on_earnings.jl")

