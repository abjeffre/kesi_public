####################
#LOAD PACKAGES #####

######################################
############ DETERMINE COMPUTER ######
using Distributed
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
@everywhere using CSV



#####################################
######## Initalize Functions ########

@everywhere include(string(pwd(), "\\functions\\utility.jl"))

######################################
#### Initalize submodules ############

@everywhere files = readdir(string(pwd(), ("\\code\\abm\\submodules")))
@everywhere for i in files  include(string(pwd(), "\\code\\abm\\submodules\\$i")) end

######################################
######### CHOOSE ABM VERSION #########

@everywhere include(string(pwd(), "\\code\\abm\\abm_cleaned.jl"))

#########################################
###########  SEASONALITY PLOT ###########

@everywhere include(string(pwd(),  "\\code\\plotting\\base_seasonality_theory_plot.jl"))

