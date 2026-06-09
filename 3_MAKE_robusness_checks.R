
##############################################################
########### ROBUSTNESS CHECK 1 - OBSERVATIONS ONLY ###########
data <- readRDS("~/kesi/data/data_kesi2025-09-22.RDS")
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
data_kesi$period_length <- 26
data_kesi$M = 2
data_kesi$N_cases = data_kesi$N
data_kesi$target_var = .1


###### RUN MODELS ########

obs_only_sector <- cmdstan_model("~/kesi/code/stan_models/obs_only_sectors.stan")
obs_only_sector_hmc <- obs_only_sector$sample(
  data = data_kesi,
  iter_sampling = 500,
  iter_warmup = 500,
  chains = 4,
  #init = pf2,
  parallel_chains = 4,
  refresh = 2
)

library(cmdstanr)
obs_only_agg <- cmdstan_model("~/kesi/code/stan_models/obs_only_aggregate.stan")

obs_only_agg_hmc <- obs_only_agg$sample(
  data = data_kesi,
  iter_sampling = 500,
  iter_warmup = 500,
  chains = 4,
  #  init = pf2,
  parallel_chains = 4,
  refresh = 2
)

#########################################
########### CALCULATRE ELASTICITY #######

post <- extract.samples(obs_only_agg_hmc)
mu_draw <- post$bgdp
a_draw  <- post$a

# --- build a FIXED (non-random) GDP baseline from posterior means ---
N <- data_kesi$N_cases
`%||%` <- function(a,b) if (!is.null(a)) a else b

gdp_true_bar   <- rowSums(data_kesi$y) 
# --- compute reference summaries ---
S_bar <- mean(log1p(gdp_true_bar))                # \overline{S}
a_bar <- mean(a_draw)
mu_bar <- mean(mu_draw)

softplus <- function(x) log1p(exp(x))
logistic <- function(x) 1/(1+exp(-x))

eta_ref <- a_bar + mu_bar * S_bar
kappa_ref <- S_bar * (logistic(eta_ref) / softplus(eta_ref))

# --- elasticity draws: linear rescale of bgdp_mu ---
gdp_mu <- kappa_ref * mu_draw

mean(gdp_mu)


##############################################################################
############# ROBUSTNESS CHECK 2 - HISTORICAL DATA ###########################

scale_gdp <-readRDS('data/gov_adjustments.RDS')
# Align Data
obs<-data$y[data$year[data$eco_ind] ==13,10]
gov<-scale_gdp[data$year[1:nrow(scale_gdp)] ==13,6]
clove_hist2<-scale_gdp[,6]*(sum(obs)/sum(gov))
data$clove_hist<- clove_hist2*data$scale_cpi[1:338]
clove_hist <- cmdstan_model("kesi/code/stan_models/clove_hist.stan")
# Sampling
hist_hmc <- clove_hist$sample(
  data = data,
  iter_sampling = 250,
  iter_warmup = 250,
  chains = 4,
  init = 0,
  parallel_chains = 4,
  refresh = 2
)

post<-extract.samples(hist_hmc)
saveRDS(post, "kesi/data/post_hist.RDS")
colnames(data$y)
dens(post$mu_gdp +post$bgdp[,6]*post$sigma_gdp, show.HPDI = .9, show.zero = T)
dens(post$mu_gdp +post$bgdp[,10]*post$sigma_gdp, show.HPDI = .9, show.zero = T)
