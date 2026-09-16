###################################################################
############ 3. ROBUSTNESS CHECKS #################################
# Reads  data/data_kesi2025-09-22.RDS, data/gov_adjustments.RDS
# Writes data/robustness_<model>.RDS, data/robustness_elasticities.csv,
#        data/post_hist.RDS, figures/robustness_clove_hist.pdf

library(rethinking)
library(cmdstanr)
library(posterior)
source("code/functions/utility.R")

CLOVE_SECTOR <- 10
CLOVE_GOV_COLUMN <- 6
ALIGNMENT_YEAR <- 13
ROBUSTNESS_PERIOD_LENGTH <- 26
ROBUSTNESS_M <- 2

fit_robustness_model <- function(stan_file, model_data) {
  model <- cmdstan_model(stan_file)
  model$sample(
    data = model_data,
    iter_sampling = SAMPLER_ITER,
    iter_warmup = SAMPLER_ITER,
    chains = SAMPLER_CHAINS,
    init = 0,
    parallel_chains = SAMPLER_CHAINS,
    max_treedepth = SAMPLER_MAX_TREEDEPTH,
    refresh = SAMPLER_REFRESH
  )
}

##############################################################
########### ROBUSTNESS CHECK 1 - OBSERVATIONS ONLY ###########
data <- readRDS("data/data_kesi2025-09-22.RDS")
if (!is.matrix(data$y)) stop("'data$y' must be a matrix (N x K).")
N <- nrow(data$y)
K <- ncol(data$y)

## helper: enforce exact length N for N-axis vectors (truncate if longer, error if shorter)
.enforce_lenN <- function(x, name) {
  if (length(x) < N) stop(sprintf("'%s' length (%d) < N from y (%d).", name, length(x), N))
  if (length(x) > N) x <- x[seq_len(N)]
  x
}
## ---- pull & align N-axis vectors ----
kesi          <- .enforce_lenN(data$kesi,          "kesi")
timber_prices <- .enforce_lenN(data$timber_prices, "timber_prices")
year          <- .enforce_lenN(data$year,          "year")
period        <- .enforce_lenN(data$period,        "period")
ramadan       <- .enforce_lenN(data$ramadan,       "ramadan")

## ---- survey version  ----
if (is.null(data$NV)) stop("'NV' (number of survey versions) missing.")
if (is.null(data$sver)) stop("'sver' vector missing.")
NV   <- as.integer(data$NV)
sver <- .enforce_lenN(data$sver, "sver")
if (any(!is.finite(sver))) stop("Non-finite values in 'sver'.")
if (any(sver < 1L | sver > NV)) stop("Some 'sver' indices are outside 1..NV.")
storage.mode(sver) <- "integer"

## ---- basic checks ----
if (any(!is.finite(kesi)))          stop("Non-finite values in 'kesi'.")
if (any(!is.finite(timber_prices))) stop("Non-finite values in 'timber_prices'.")
if (any(year   < 1 | year   > data$L))  stop("'year' outside [1, L].")
if (any(period < 1 | period > data$LP)) stop("'period' outside [1, LP].")
storage.mode(year)    <- "integer"
storage.mode(period)  <- "integer"
storage.mode(ramadan) <- "integer"

## ---- start payload (y untouched; N comes from y) ----
data_kesi <- list(
  N  = as.integer(N),
  K  = as.integer(K),
  L  = as.integer(data$L),
  LP = as.integer(data$LP),
  kesi          = as.integer(kesi),
  y             = data$y,                    # DO NOT touch y
  timber_prices = as.numeric(timber_prices),
  year          = as.integer(year),
  period        = as.integer(period),
  ramadan       = as.integer(ramadan),
  ## survey version
  NV   = NV,
  sver = sver,
  ## time kernels (unchanged shapes)
  DmatX         = data$DmatX,               # [L,L]
  P1            = data$P1,                  # length LP
  period_length = data$period_length,
  ## link stability
  softplus_alpha = data$softplus_alpha
)

## ---- direct env GPs  ----
if (is.null(data$DM)) stop("'data$DM' missing.")
data_kesi$DM <- as.integer(data$DM)
if (data_kesi$DM == 0L) {
  data_kesi$NP       <- 0L
  data_kesi$DmatDEnv <- if (!is.null(data$DmatDEnv)) data$DmatDEnv else array(0, dim = c(0,0,0))
  data_kesi$denv_ind <- matrix(0L, nrow = 0L, ncol = N)
} else {
  if (is.null(data$NP))        stop("'data$NP' missing while DM>0.")
  if (is.null(data$DmatDEnv))  stop("'data$DmatDEnv' missing while DM>0.")
  if (is.null(data$denv_ind))  stop("'data$denv_ind' missing while DM>0.")
  
  data_kesi$NP       <- as.integer(data$NP)
  data_kesi$DmatDEnv <- data$DmatDEnv   # KEEP exact 3D shape; no permute
  ## denv_ind must be [DM x N]; truncate columns if longer; error if shorter
  if (is.matrix(data$denv_ind)) {
    if (nrow(data$denv_ind) != data_kesi$DM) stop("nrow(denv_ind) != DM.")
    if (ncol(data$denv_ind) < N)
      stop(sprintf("denv_ind has %d cols < N (%d).", ncol(data$denv_ind), N))
    denv_ind_mat <- data$denv_ind[, seq_len(N), drop = FALSE]
  } else if (is.list(data$denv_ind)) {
    if (length(data$denv_ind) != data_kesi$DM) stop("length(denv_ind) != DM (list).")
    denv_ind_mat <- matrix(NA_integer_, nrow = data_kesi$DM, ncol = N)
    for (m in seq_len(data_kesi$DM)) {
      if (length(data$denv_ind[[m]]) < N)
        stop(sprintf("denv_ind[[%d]] shorter than N.", m))
      denv_ind_mat[m, ] <- as.integer(data$denv_ind[[m]][seq_len(N)])
    }
  } else {
    stop("denv_ind must be a matrix [DM x N_full] or list(DM) of length-N_full vectors.")
  }
  
  storage.mode(denv_ind_mat) <- "integer"
  if (any(!is.finite(denv_ind_mat))) stop("denv_ind has non-finite values in first N columns.")
  if (any(denv_ind_mat < 1L | denv_ind_mat > data_kesi$NP))
    stop("denv_ind indices outside 1..NP in first N columns.")
  
  data_kesi$denv_ind <- denv_ind_mat
}
data_kesi$period_length <- ROBUSTNESS_PERIOD_LENGTH
data_kesi$M <- ROBUSTNESS_M
data_kesi$N_cases <- data_kesi$N

###### RUN MODELS ########

obs_only_sector_hmc <- fit_robustness_model("code/stan_models/obs_only_sectors.stan", data_kesi)
saveRDS(extract.samples(obs_only_sector_hmc), "data/robustness_obs_only_sectors.RDS")

obs_only_agg_hmc <- fit_robustness_model("code/stan_models/obs_only_aggregate.stan", data_kesi)
saveRDS(extract.samples(obs_only_agg_hmc), "data/robustness_obs_only_aggregate.RDS")

# Full-model variants: same data as the main model, older name for the CPI scale
data_variants <- data
data_variants$scale_gdp <- data$scale_cpi

model_sectors_hmc <- fit_robustness_model("code/stan_models/model_sectors.stan", data_variants)
saveRDS(extract.samples(model_sectors_hmc), "data/robustness_model_sectors.RDS")

model_aggregate_hmc <- fit_robustness_model("code/stan_models/model_aggregate.stan", data_variants)
saveRDS(extract.samples(model_aggregate_hmc), "data/robustness_model_aggregate.RDS")

#########################################
########### CALCULATE ELASTICITY ########

post <- extract.samples(obs_only_agg_hmc)
mu_draw <- post$bgdp
a_draw  <- post$a

# --- build a FIXED (non-random) GDP baseline from posterior means ---
gdp_true_bar   <- rowSums(data_kesi$y)
S_bar <- mean(log1p(gdp_true_bar))
a_bar <- mean(a_draw)
mu_bar <- mean(mu_draw)

softplus <- function(x) log1p(exp(x))
logistic <- function(x) 1/(1+exp(-x))

eta_ref <- a_bar + mu_bar * S_bar
kappa_ref <- S_bar * (logistic(eta_ref) / softplus(eta_ref))

# --- elasticity draws: linear rescale of bgdp_mu ---
gdp_mu <- kappa_ref * mu_draw

write.csv(data.frame(model = "obs_only_aggregate",
                     elasticity_mean = mean(gdp_mu),
                     elasticity_lower_90 = PI(gdp_mu, .9)[1],
                     elasticity_upper_90 = PI(gdp_mu, .9)[2]),
          "data/robustness_elasticities.csv", row.names = FALSE)

##############################################################################
############# ROBUSTNESS CHECK 2 - HISTORICAL DATA ###########################

scale_gdp <- readRDS("data/gov_adjustments.RDS")
# Align the official series to the survey earnings in the overlap year
obs <- data$y[data$year[data$eco_ind] == ALIGNMENT_YEAR, CLOVE_SECTOR]
gov <- scale_gdp[data$year[1:nrow(scale_gdp)] == ALIGNMENT_YEAR, CLOVE_GOV_COLUMN]
clove_hist2 <- scale_gdp[, CLOVE_GOV_COLUMN] * (sum(obs) / sum(gov))
# The official series covers the first nrow(scale_gdp) periods of the timeline;
# only the imputed periods (est_ind) are replaced, observed survey periods stay.
clove_hist_timeline <- clove_hist2 * data$scale_cpi[1:nrow(scale_gdp)]
data$clove_hist <- clove_hist_timeline[data$est_ind]
data$clove_sector <- CLOVE_SECTOR

hist_hmc <- fit_robustness_model("code/stan_models/clove_hist.stan", data)

post <- extract.samples(hist_hmc)
saveRDS(post, "data/post_hist.RDS")

pdf("figures/robustness_clove_hist.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
dens(post$mu_gdp + post$bgdp[, 6] * post$sigma_gdp, show.HPDI = .9, show.zero = TRUE,
     xlab = paste("Effect of", colnames(data$y)[6]))
dens(post$mu_gdp + post$bgdp[, CLOVE_SECTOR] * post$sigma_gdp, show.HPDI = .9, show.zero = TRUE,
     xlab = paste("Effect of", colnames(data$y)[CLOVE_SECTOR]))
dev.off()
