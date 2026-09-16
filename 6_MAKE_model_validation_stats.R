###################################################################
############ 6. MODEL VALIDATION - RECOVERY FITS ##################
# Reads  data/sweeps/*/abm/*.csv and sweep_list.csv (from 5_MAKE)
# Writes data/sweeps/*/stan/*.csv, figures/prediction_error.pdf,
#        figures/combined_plot.pdf
# KESI_CORES caps the parallel workers; KESI_SMOKE fits the first few cells only.

library(readr)
library(rethinking)
library(posterior)
library(cmdstanr)
library(ggplot2)
library(cowplot)
source("code/functions/utility.R")

source("code/sweeps/weather_on_earnings_recover.R")
source("code/sweeps/earnings_on_kesi_recover.R")

source("code/plotting/weather_on_earnings_sweep_heatmap.R")
source("code/plotting/earnings_on_kesi_sweep_heatmap.R")
