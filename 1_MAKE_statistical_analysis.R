###################################################################
############ 1. FIT THE MAIN MODEL ################################
# Reads  data/data_kesi2025-09-22.RDS
# Writes data/full_hmc.RDS (posterior draws, ~1.5 GB at full settings)

library(rethinking)
library(cmdstanr)
library(posterior)
source("code/functions/utility.R")

data <- readRDS("data/data_kesi2025-09-22.RDS")

full_24_model <- cmdstan_model("code/stan_models/main_model.stan")

full_hmc <- full_24_model$sample(
  data = data,
  iter_sampling = SAMPLER_ITER,
  iter_warmup = SAMPLER_ITER,
  chains = SAMPLER_CHAINS,
  init = 0,
  parallel_chains = SAMPLER_CHAINS,
  max_treedepth = SAMPLER_MAX_TREEDEPTH,
  refresh = SAMPLER_REFRESH
)
post <- extract.samples(full_hmc)
saveRDS(post, "data/full_hmc.RDS")
