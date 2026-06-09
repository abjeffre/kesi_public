####################################################################
################## MAKE THE BASE PLOT ##############################

dens2 <-function (x, adj = 0.5, norm.comp = FALSE, main = "", show.HPDI = FALSE, 
                  show.zero = FALSE, rm.na = TRUE, add = FALSE, col_fill = "grey", ...) 
{
  if (inherits(x, "data.frame")) {
    n <- ncol(x)
    cnames <- colnames(x)
    set_nice_margins()
    par(mfrow = make.grid(n))
    for (i in 1:n) {
      dens(x[, i], adj = adj, norm.comp = norm.comp, show.HPDI = show.HPDI, 
           show.zero = TRUE, xlab = cnames[i], ...)
    }
  }
  else {
    if (rm.na == TRUE) 
      x <- x[!is.na(x)]
    thed <- density(x, adjust = adj)
    if (add == FALSE) {
      set_nice_margins()
      plot(thed, main = main, ...)
    }
    else lines(thed$x, thed$y, ...)
    if (show.HPDI != FALSE) {
      hpd <- HPDI(x, prob = show.HPDI)
      shade(thed, hpd, col = col_fill)
    }
    if (norm.comp == TRUE) {
      mu <- mean(x)
      sigma <- sd(x)
      curve(dnorm(x, mu, sigma), col = "white", lwd = 2, 
            add = TRUE)
      curve(dnorm(x, mu, sigma), add = TRUE)
    }
    if (show.zero == TRUE) {
      lines(c(0, 0), c(0, max(thed$y) * 2), lty = 2)
    }
  }
}

dotchart2<-function (x, labels = NULL, groups = NULL, gdata = NULL, offset = 1/8, 
                     ann = par("ann"), xaxt = par("xaxt"), frame.plot = TRUE, 
                     log = "", cex = par("cex"), pt.cex = cex, pch = 21, gpch = 21, 
                     bg = par("bg"), color = par("fg"), gcolor = par("fg"), lcolor = "gray", 
                     xlim = range(x[is.finite(x)]), main = NULL, xlab = NULL, tcolor = "black",
                     ylab = NULL, ...) 
{
  opar <- par("mai", "mar", "mgp", "cex", "yaxs")
  on.exit(par(opar))
  par(cex = cex, yaxs = "i")
  if (!is.numeric(x)) 
    stop("'x' must be a numeric vector or matrix")
  n <- length(x)
  if (is.matrix(x)) {
    if (is.null(labels)) 
      labels <- rownames(x)
    if (is.null(labels)) 
      labels <- as.character(seq_len(nrow(x)))
    labels <- rep_len(labels, n)
    if (is.null(groups)) 
      groups <- col(x, as.factor = TRUE)
    glabels <- levels(groups)
  }
  else {
    if (is.null(labels)) 
      labels <- names(x)
    glabels <- if (!is.null(groups)) 
      levels(groups)
    if (!is.vector(x)) {
      warning("'x' is neither a vector nor a matrix: using as.numeric(x)")
      x <- as.numeric(x)
    }
  }
  plot.new()
  linch <- if (!is.null(labels)) 
    max(strwidth(labels, "inch"), na.rm = TRUE)
  else 0
  if (is.null(glabels)) {
    ginch <- 0
    goffset <- 0
  }
  else {
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- offset
  }
  nmai <- opar[["mai"]]
  if (ann) 
    nm.2 <- nmai[2L]
  if (!(is.null(labels) && is.null(glabels))) {
    yi <- if (is.null(ylab) || !ann) 
      0
    else offset
    nm.2 <- nmai[4L] + max(yi + linch + goffset, ginch) + 
      1/16
    if (nmai[2L] < nm.2) {
      nmai[2L] <- nm.2
      par(mai = nmai)
    }
  }
  if (is.null(groups)) {
    o <- seq_len(n)
    y <- o
    ylim <- c(0, n + 1)
  }
  else {
    o <- sort.list(as.numeric(groups), decreasing = TRUE)
    x <- x[o]
    groups <- groups[o]
    color <- rep_len(color, length(groups))[o]
    lcolor <- rep_len(lcolor, length(groups))[o]
    pch <- rep_len(pch, length(groups))[o]
    of.1 <- cumsum(c(0, diff(as.numeric(groups)) != 0))
    y <- seq_len(n) + 2 * of.1
    ylim <- range(0, y + 2)
  }
  plot.window(xlim = xlim, ylim = ylim, log = log)
  lheight <- par("csi")
  if (!is.null(labels)) {
    loffset <- (linch + 0.1)/lheight
    mtext(labels[o], side = 2, line = loffset, at = y, adj = 0, 
          col = tcolor, las = 2, cex = cex, ...)
  }
  abline(h = y, lty = "dotted", col = lcolor)
  points(x, y, pch = pch, col = color, bg = bg, cex = pt.cex/cex)
  if (!is.null(groups)) {
    gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
                         2) - 1)
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- (max(linch + offset, ginch, na.rm = TRUE) + 
                  1/16)/lheight
    mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
          col = gcolor, las = 2, cex = cex, ...)
    if (!is.null(gdata)) {
      abline(h = gpos, lty = "dotted")
      points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
             cex = pt.cex/cex, ...)
    }
  }
  axis(1, xaxt = xaxt)
  if (frame.plot) 
    box()
  if (ann) {
    title(main = main, xlab = xlab, ...)
    mgp <- par("mgp")
    par(mgp = c(max(mgp[1], nm.2/lheight - 1.5), mgp[-1]))
    title(ylab = ylab, ...)
  }
  invisible()
}

##### LOAD DATA ############
post <- readRDS('data/full_hmc.RDS')
# --- helper: pick param name ---
mu_draw <- if ("bgdp_mu" %in% names(post)) post$bgdp_mu else post$mu_gdp
a_draw  <- post$a

# --- build a FIXED (non-random) GDP baseline from posterior means ---
N_env <- data$N_env; K <- data$K
est_ind <- data$est_ind; eco_ind <- data$eco_ind
`%||%` <- function(a,b) if (!is.null(a)) a else b
clove_hist <- data$clove_hist %||% NULL

gdp_impute_bar <- apply(post$gdp_impute, c(2,3), mean)  # [N_est,K]
gdp_true_bar   <- apply(post$gdp_true,   c(2,3), mean)  # [N_eco,K]

gm_fixed <- matrix(NA_real_, N_env, K)
gm_fixed[est_ind, ] <- gdp_impute_bar
gm_fixed[eco_ind, ] <- gdp_true_bar
if (!is.null(clove_hist) && K == 10) gm_fixed[est_ind, 10] <- clove_hist

# --- compute reference summaries ---
S_bar <- mean(rowSums(log1p(gm_fixed)))                # \overline{S}
C_bar <- mean(rowSums(gm_fixed/(1+gm_fixed)))          # \overline{C}
a_bar <- mean(a_draw)
mu_bar <- mean(mu_draw)

softplus <- function(x) log1p(exp(x))
logistic <- function(x) 1/(1+exp(-x))

eta_ref <- a_bar + mu_bar * S_bar
kappa_ref <- C_bar * (logistic(eta_ref) / softplus(eta_ref))

# --- elasticity draws: linear rescale of bgdp_mu ---
gdp_mu <- kappa_ref * mu_draw


###### MAKE PLOTS ######
#post <- post
pdf(file ="figures/base_relationships.pdf", width = 16, height = 6)

par(mfrow = c(1, 3), mar = c(6,4.2,2, 1))
dens2(gdp_mu, show.HPDI = .90,
      xlab = "Elasticity of Earnings on Illegal Activity",
      col_fill = "#7fc5bd", col = "black", show.zero = T, cex.lab  = 1.7, cex.axis = 1.5, cex = 1.3 )
title("A.", adj = 0)

if(!"gdp_merged" %in% names(post)) post$gdp_merged <-  abind::abind(post$gdp_impute, post$gdp_true, along = 2)
gdp_impute<-colMeans(apply(post$gdp_merged, 2, rowSums))
seq = seq(60,140, length.out = 100)
sim<- function(seq){
  gdp<-matrix(NA, ncol = data$K, nrow=nrow(post$mu_gdp))
  for(k in 1:data$K){
    gdp[,k] <- (post$mu_gdp+post$sigma_gdp*post$bgdp[,k])*log1p(seq * mean(colMeans(post$gdp_merged[,,k])/gdp_impute))
  }
  rowSums(gdp)
  rnbinom(1000, mu = with(post, 
                          log1p(exp(a +
                                      rowSums(gdp) +
                                      rowMedians(kDEnv[,2,]) +
                                      rowMedians(kDEnv[,1,]) +
                                      rowMedians(kp[,]) +
                                      rowMedians(kX[,]) +
                                      rowMedians(br_k[,]) +
                                      btimber*mean(data$timber_prices)))), post$omega)
}

plot(data$kesi ~ gdp_impute[1:length(data$kesi)] , xlab = "Household Earnings (Two week period) ", ylab = "Illegal Activity",  pch = 16,
     col =col.alpha("#606161", .5), xlim = c(63, 137), cex.lab = 1.6, cex = 1.7, cex.axis = 1.5)
title("B.", adj = 0)
ests <- sapply(seq, sim)
rethinking::shade(apply(ests, 2, PI, .68), seq, col = col.alpha("#7fc5bd", .2))
lines(colMeans(ests) ~ (seq), lwd = 4, col  = "#006c66")

lines(c(0, 160), c(colMeans(ests)[1], colMeans(ests)[1]), lwd = 2, col  = "#6C0006", lty = 2)
lines(c(0, 160), c(colMeans(ests)[100], colMeans(ests)[100]), lwd = 2, col  = "#6C0006", lty = 2)
#points(data$kesi ~ (colMeans(post$gdp_merged)), pch = 16,  col = col.alpha("#606161", .1))
text(x = 132, y = c(colMeans(ests)[1])+.4, col  = "#6C0006", cex = 1.5, # Coordinates
     label = paste0("~ ", round(colMeans(ests)[1],2)))
text(x = 132, y = c(colMeans(ests)[100])-.4,  col  = "#6C0006", cex = 1.5, # Coordinates
     label = paste0("~ ", round(colMeans(ests)[100],2)))

text(x = 120, y = 2.40,  col  = "#6C0006", # Coordinates
     label = paste0("~ ", round((1-(colMeans(ests)[100]/colMeans(ests)[1]))*100,0) ,"% Decrease"), cex = 1.5)
# 
# library(shape)
# shape::Arrows(
#   132, colMeans(ests)[100] + 0.3,   # start (x0, y0) — above
#   132, colMeans(ests)[1]  - 0.3,    # end   (x1, y1) — below
#   code = 2,
#   arr.type = "triangle",
#   arr.length = 0.25,
#   arr.width  = 0.25,
#   col = "#6C0006"
# )
# 

library(shape)
shape::Arrows(
  132, colMeans(ests)[1] - 0.3,   # start (x0, y0) — above
  132, colMeans(ests)[100]  + 0.3,    # end   (x1, y1) — below
  code = 2,
  arr.type = "triangle",
  arr.length = 0.25,
  arr.width  = 0.25,
  col = "#6C0006"
)



#dev.off()
#pdf(file ="figures/catapiller.pdf", width = 7, height = 6)
# 
# df <- NULL
# for(i in 1:data$K){
#   df<-cbind(df, post$bgdp[,i]*post$sigma_gdp + post$mu_gdp)
# }

# === Elasticities from the original (softplus) model ===
# Works with: post <- rethinking::extract.samples(fit)
# and the same data list used for Stan. Optional: data_list$clove_hist (length N_est)
elasticities_from_post <- function(
    post, data_list,
    draws = NULL,
    return_obs_level = FALSE,
    marginalize = FALSE,     # if TRUE: ignore kp, kX, direct env, ramadan, timber (GDP-only eta)
    period_length = 26.0     # your seasonal period (matches Stan)
) {
  # ---- quick aliases ----
  N_cases <- data_list$N_cases
  N_est   <- data_list$N_est
  N_eco   <- data_list$N_eco
  N_env   <- data_list$N_env
  K       <- data_list$K
  DM      <- data_list$DM
  LP      <- data_list$LP
  L       <- data_list$L
  
  period  <- data_list$period    # length N_env
  year    <- data_list$year      # length N_env
  ramadan <- data_list$ramadan   # length N_env in {1,2}
  denv_ind<- data_list$denv_ind  # DM x N_env (indices 1..NP_m)
  eco_ind <- data_list$eco_ind   # length N_eco (env row per eco obs)
  est_ind <- data_list$est_ind   # length N_est (env row per "est" slot)
  timber  <- data_list$timber_prices
  DmatX   <- data_list$DmatX
  DmatDEnv<- data_list$DmatDEnv
  P1      <- data_list$P1
  
  clove_hist <- data_list$clove_hist %||% NULL  # optional N_est vector
  
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  as_mat <- function(x) { # safe coercion
    x <- as.matrix(x); storage.mode(x) <- "double"; x
  }
  
  # ---- robust matrix getter for DmatDEnv ----
  # Uses the length of zDEnv[s,m,] (NP_m) to slice a square NP_m x NP_m matrix.
  get_mat_m <- function(D, m, NP_m) {
    dx <- dim(D)
    if (length(dx) == 3L) {
      stopifnot(dx[1] >= NP_m, dx[2] >= NP_m)
      return(as_mat(D[1:NP_m, 1:NP_m, m]))
    } else {
      M <- as_mat(D[[m]])
      stopifnot(nrow(M) >= NP_m, ncol(M) >= NP_m)
      return(as_mat(M[1:NP_m, 1:NP_m]))
    }
  }
  
  # ---- kernels (mirror Stan) ----
  cov_GPL2 <- function(D, sq_alpha, sq_rho, delta) {
    Dm <- as_mat(D)
    K  <- sq_alpha * exp(-sq_rho * (Dm^2))
    diag(K) <- diag(K) + delta
    K
  }
  gp_periodic_cov_R <- function(P, alpha, rho, period) {
    LP <- length(P)
    K  <- matrix(0, LP, LP)
    for (i in 1:LP) {
      for (j in i:LP) {
        arg <- pi * (P[i] - P[j]) / period
        val <- (alpha^2) * exp(-2 * (sin(arg)^2) / (rho^2))
        K[i,j] <- val; K[j,i] <- val
      }
    }
    diag(K) <- diag(K) + 1e-9
    K
  }
  
  softplus <- function(x) log1p(exp(x))   # alpha=1
  logistic <- function(x) 1/(1+exp(-x))
  
  # ---- choose draws ----
  S_all <- length(post$a)
  if (is.null(draws)) draws <- 1:S_all
  
  # ---- builders (per draw) ----
  build_gdp_merged <- function(s) {
    gm <- matrix(NA_real_, N_env, K)
    gm[est_ind, ] <- post$gdp_impute[s,,]  # N_est x K
    gm[eco_ind, ] <- post$gdp_true[s,,]    # N_eco x K
    # Optional sector-10 overwrite for imputed rows:
    if (!is.null(clove_hist) && K == 10) {
      gm[est_ind, 10] <- clove_hist
    }
    gm
  }
  
  build_kp <- function(s) {
    Kp <- gp_periodic_cov_R(P1, post$scalesq[s], post$length_scale[s], period_length)
    drop(chol(Kp) %*% post$kp_hat[s,])
  }
  build_kX <- function(s) {
    KX <- cov_GPL2(DmatX, post$etasq[s], post$rhosq[s], post$sigmaX[s])
    drop(chol(KX) %*% post$zX[s,])
  }
  
  build_direct_env_effect <- function(s) {
    out <- numeric(N_env)
    for (m in 1:DM) {
      z_m  <- post$zDEnv[s, m, ]        # length NP_m
      NP_m <- length(z_m)
      Dm   <- get_mat_m(DmatDEnv, m, NP_m)
      K_m  <- cov_GPL2(Dm, post$etasq_DEnv[s,m], post$rhosq_DEnv[s,m], post$sigma_DEnv[s,m])
      kvec <- drop(chol(K_m) %*% z_m)   # NP_m
      out  <- out + kvec[ denv_ind[m,] ]
    }
    out
  }
  
  # ---- allocate outputs ----
  if (return_obs_level) {
    E_arr <- array(NA_real_, dim = c(length(draws), N_cases, K))
  } else {
    E_mat <- matrix(NA_real_, nrow = length(draws), ncol = K)
  }
  
  # ---- main loop ----
  for (idx in seq_along(draws)) {
    s <- draws[idx]
    
    # effective GDP slope vector (matches your first model)
    beta_eff <- post$bgdp[s,] * post$sigma_gdp[s] + post$mu_gdp[s]
    
    gdp_merged <- build_gdp_merged(s)       # N_env x K
    log1p_g    <- log1p(gdp_merged)         # N_env x K
    dot_g      <- as.numeric(log1p_g %*% beta_eff)  # N_env
    
    # pieces for eta
    a    <- post$a[s]
    bt   <- post$btimber[s]
    kp_i <- if (marginalize) rep(0, N_env) else build_kp(s)[period]
    kX_i <- if (marginalize) rep(0, N_env) else build_kX(s)[year]
    envd <- if (marginalize) rep(0, N_env) else build_direct_env_effect(s)
    brk  <- if (marginalize) c(0,0) else post$br_k[s,]
    
    eta <- a + kp_i + dot_g + kX_i + brk[ramadan] + envd + bt * timber
    
    lambda  <- softplus(eta)            # α = 1
    sig_eta <- logistic(eta)
    frac    <- sig_eta / lambda         # length N_env
    
    ratio   <- gdp_merged / (1 + gdp_merged)  # N_env x K
    E_obsk  <- sweep(ratio, 1, frac, `*`)
    E_obsk  <- sweep(E_obsk, 2, beta_eff, `*`)
    
    if (return_obs_level) {
      E_arr[idx,,] <- E_obsk
    } else {
      E_mat[idx,] <- colMeans(E_obsk)
    }
  }
  
  if (return_obs_level) {
    ame <- apply(E_arr, c(1,3), mean)  # [draw, k]
    list(
      elasticities_draw_i_k = E_arr,
      AME_draw_k = ame,
      AME_summary = data.frame(
        sector = 1:K,
        mean = colMeans(ame),
        sd   = apply(ame, 2, sd),
        lo95 = apply(ame, 2, quantile, 0.025),
        hi95 = apply(ame, 2, quantile, 0.975)
      )
    )
  } else {
    list(
      AME_draw_k = E_mat,
      AME_summary = data.frame(
        sector = 1:K,
        mean = colMeans(E_mat),
        sd   = apply(E_mat, 2, sd),
        lo95 = apply(E_mat, 2, quantile, 0.025),
        hi95 = apply(E_mat, 2, quantile, 0.975)
      )
    )
  }
}




# 2) Marginalize everything but GDP (safer if your distance arrays are weird):
res_m <- elasticities_from_post(post, data, return_obs_level = FALSE, marginalize = TRUE)
df<-res_m$AME_draw_k


# STEP 1 GET MEANS

set.seed(123)


mu  = rev(colMeans(df))

names(mu) <- rev(colnames(data$y))
names(mu) <- c("Other Services", "Cloves", "Trade and Transport", "Government",
               "Seaweed", "Fishing", "Manufacturing", "Forestry", "Construction",
               "Livestock", "Farming")

# Calculate error bars for each column
errors <- apply(df, 2, PI, .89)[, ncol(df):1]

# Create dot chart
dotchart2(mu,
          cex = 1.0,
          xlab = "Elasticity of Earnings on Illegal Activity",
          pch = 16,
          xlim = c(-.5, .5), 
          col  = "#006c66",
          tcolor = "black",
          cex.lab = 1.25)

# Add error bars
for (i in 1:ncol(df)) {
  #points(x = errors[,i ],
  #     y = rep(i, 2),
  #     pch = 19,
  #     col = i,
  #     xlim = range(df) + c(-1, 1) * max(errors),
  #     ylim = c(0.5, ncol(df) + 0.5))
  segments(x0 = errors[1,i],
           x1 = errors[2,i],
           y0 = rep(i, 1),
           y1 = rep(i, 1),
           lwd = 2,
           col  = "#006c66")
}
abline(v = 0, col = "black", lty = 2)
title("C.", adj = 0)
dev.off()


#################################################
################ OND ############################
## =========================
## 3-PANEL OND FIGURE (robust full script)
## =========================

post <- readRDS("data/full_hmc.RDS")
data <- readRDS("data/data_kesi2025-09-22.RDS")

## -------- palette ----------
col_teal   <- "#1B9088"
col_teal_d <- "#0D5E59"
col_light  <- "#E8F4F3"
col_border <- "#BCCAC8"

## -------- utilities ----------
softplus <- function(x, alpha = 1) log1p(exp(alpha * x)) / alpha
hpdi <- function(x, prob = 0.90) {
  x <- sort(as.numeric(x)); n <- length(x); m <- max(1, ceiling(prob*n))
  j <- which.min(x[m:n] - x[1:(n-m+1)])
  c(x[j], x[j+m-1])
}
.safe_range <- function(x) {
  xf <- x[is.finite(x)]
  if (length(xf)) range(xf) else c(-1, 1)
}
.pretty_at <- function(x) pretty(x[is.finite(x)], n = 5)

gp_periodic_cov_R <- function(P, alpha, rho, period_len = 26.0) {
  LP <- length(P); Kmat <- matrix(0, LP, LP)
  for (i in 1:LP) for (j in i:LP) {
    a <- pi*(P[i]-P[j])/period_len
    v <- (alpha^2) * exp(-2 * (sin(a)^2) / (rho^2))
    Kmat[i,j] <- v; Kmat[j,i] <- v
  }
  diag(Kmat) <- diag(Kmat) + 1e-9
  Kmat
}

## -------- lag window helpers (OND-aware) ----------
month_from_period <- function(p, LP = 26L) floor(12*((as.numeric(p)-0.5)/LP)) + 1L

# Return the set of months impacted by an OND shift for a given lag spec.
# lag_spec can be "inst" (offsets 0:0), or an integer L meaning offsets (L-1):(L).
impacted_months_for_lag <- function(lag_spec, ond_months = c(10L,11L,12L)) {
  add_mod <- function(m, off) ((m - 1L + off) %% 12L) + 1L
  if (identical(lag_spec, "inst")) {
    return(sort(unique(ond_months)))
  } else {
    L <- as.integer(lag_spec)
    offs <- (L - 1L):L
    res <- unlist(lapply(ond_months, function(m0) add_mod(m0, offs)))
    return(sort(unique(res)))
  }
}

# Build per-variable row-index masks (over ENV rows used in estimation) that should be shifted
# when OND is perturbed. Also returns the union of all impacted months for summary windows.
make_lag_row_masks <- function(per, est_idx, target_vars, rn) {
  # map from env var name -> months impacted
  # recognized patterns: "ind_moisture" (instant), "ind_lag{L}_rain" where L in 2..6
  var_months <- list()
  for (v in target_vars) {
    if (v == "ind_moisture") {
      var_months[[v]] <- impacted_months_for_lag("inst")
    } else {
      m <- regexec("^ind_lag([0-9]+)_rain$", v)
      reg <- regmatches(v, m)[[1]]
      if (length(reg)) {
        L <- as.integer(reg[2])
        var_months[[v]] <- impacted_months_for_lag(L)
      } else {
        # fallback: treat unknown as contemporaneous
        var_months[[v]] <- impacted_months_for_lag("inst")
      }
    }
  }
  months_all <- sort(unique(unlist(var_months)))
  
  # convert months -> row indices (restrict to estimation rows so we only shift those)
  per_est_month <- month_from_period(per[est_idx])
  masks <- lapply(var_months, function(mm) est_idx[ per_est_month %in% mm ])
  
  list(masks = masks, months_all = months_all)
}

## -------- sector-level draws (lag-aware shifting & impact window) ----------
build_sector_bin_draws <- function(post, data,
                                   bin_seq = -4:4,
                                   target_vars = c("ind_moisture",
                                                   "ind_lag2_rain","ind_lag3_rain",
                                                   "ind_lag4_rain","ind_lag5_rain","ind_lag6_rain"),
                                   draws = 200,
                                   env_mode = c("mean","full")) {
  env_mode <- match.arg(env_mode)
  stopifnot(data$LP == 26L)
  
  N_env   <- data$N_env
  N_est   <- data$N_est
  N_eco   <- data$N_eco
  K       <- data$K
  M       <- nrow(data$env_ind)
  NP      <- dim(data$DmatEnv)[2]
  est_idx <- data$est_ind
  eco_idx <- data$eco_ind
  per     <- as.integer(data$period)   # 1..26
  ram     <- as.integer(data$ramadan)  # 1..2
  scg     <- as.numeric(data$scale_gdp); scg[!is.finite(scg) | scg <= 0] <- 1
  P1      <- as.numeric(data$P1)
  rn      <- rownames(data$env_ind)
  
  # lag-aware: which rows to shift per var; what months to summarize over
  lag_info <- make_lag_row_masks(per, est_idx, target_vars, rn)
  var_masks <- lag_info$masks
  impact_months <- lag_info$months_all
  impact_rows_env <- which(month_from_period(per) %in% impact_months)
  impact_mask_env <- seq_len(N_env) %in% impact_rows_env
  
  sec_names <- colnames(data$y); if (is.null(sec_names)) sec_names <- paste0("sector_", seq_len(K))
  
  S_all <- length(post$a)
  use   <- seq_len(min(draws, S_all))
  
  arr <- array(NA_real_, dim = c(length(bin_seq), K, length(use)),
               dimnames = list(bin = as.character(bin_seq), sector = sec_names, draw = NULL))
  
  # ---- robust extractors ----
  get_mu_EnvP_draw <- function(post, s, M, K, NP) {
    x <- post$mu_EnvP; dx <- dim(x)
    if (length(dx) == 4L) array(x[s,,, , drop = FALSE], dim = c(M, K, NP))
    else if (length(dx) == 3L) array(x, dim = c(M, K, NP))
    else stop("mu_EnvP dims: ", paste(dx, collapse="x"))
  }
  get_LEnv_draw <- function(post, s, M, K, NP) {
    x <- post$LEnv_sigmaP; dx <- dim(x)
    if (length(dx) == 5L) array(x[s,,,, , drop=FALSE], dim = c(M,K,NP,NP))
    else if (length(dx) == 4L) array(x, dim = c(M,K,NP,NP))
    else stop("LEnv_sigmaP dims: ", paste(dx, collapse="x"))
  }
  get_cov_draw <- function(post, s, M, K) {
    x <- post$cov_EnvP; dx <- dim(x)
    if (length(dx) == 3L) array(x[s,, , drop=FALSE], dim = c(M,K))
    else if (length(dx) == 2L) array(x, dim = c(M,K))
    else stop("cov_EnvP dims: ", paste(dx, collapse="x"))
  }
  get_zEnv_draw <- function(post, s, P, M, K, NP) {
    x <- post$zEnvP; dx <- dim(x)
    if (length(dx) == 5L) array(x[s,,,, , drop=FALSE], dim = c(P,M,K,NP))
    else if (length(dx) == 4L) array(x, dim = c(P,M,K,NP))
    else stop("zEnvP dims: ", paste(dx, collapse="x"))
  }
  get_br_draw <- function(post, s, K) {
    x <- post$br; dx <- dim(x)
    if (length(dx) == 3L) array(x[s,, , drop = FALSE], dim = c(K,2))
    else if (length(dx) == 2L && dx[1]==K && dx[2]==2) array(x, dim = c(K,2))
    else stop("br dims: ", paste(dx, collapse="x"))
  }
  get_sigma_r_draw <- function(post, s, K) {
    x <- post$sigma_r; dx <- dim(x)
    if (length(dx) == 2L) as.numeric(x[s, ])
    else if (length(x) == K) as.numeric(x)
    else stop("sigma_r dims unexpected.")
  }
  get_gdp_true_draw <- function(post, s, N_eco, K) {
    x <- post$gdp_true; dx <- dim(x)
    if (length(dx) == 3L) array(x[s,, , drop = FALSE], dim = c(N_eco, K))
    else if (length(dx) == 2L && dx[1]==N_eco && dx[2]==K) array(x, dim = c(N_eco, K))
    else stop("gdp_true dims: ", paste(dx, collapse="x"))
  }
  
  for (d in seq_along(use)) {
    s <- use[d]
    
    # seasonal sector kernel
    kbp <- matrix(NA_real_, K, data$LP)
    for (k in 1:K) {
      Kper <- gp_periodic_cov_R(P1, post$scalesq_b[s,k], post$length_scale_b[s,k], 26.0)
      kbp[k, ] <- drop(chol(Kper) %*% post$kbp_hat[s, k, ])
    }
    
    mu_EnvP <- get_mu_EnvP_draw(post, s, M, K, NP)
    
    # --- FULL EnvP reconstruction bits (only if requested) ---
    if (env_mode == "full") {
      LEnv   <- get_LEnv_draw(post, s, M, K, NP)     # [M,K,NP,NP]
      covEP  <- get_cov_draw(post, s, M, K)          # [M,K]
      zEnv   <- get_zEnv_draw(post, s, data$LP, M, K, NP)  # [LP,M,K,NP]
      # helper: value for (p,m,k,j)
      kEnvP_val <- function(p_i, m, k, j) {
        v <- zEnv[p_i, m, k, ] * covEP[m, k] + mu_EnvP[m, k, ]
        out <- drop( LEnv[m, k, , ] %*% v )
        out[j]
      }
    }
    
    br       <- get_br_draw(post, s, K)
    sigma_r  <- get_sigma_r_draw(post, s, K)
    gdp_true <- get_gdp_true_draw(post, s, N_eco, K)
    
    # ---- impute from env_ind (mean vs full) ----
    impute_from_env <- function(env_ind_mat) {
      inner <- matrix(0, N_est, K)
      for (i in 1:N_est) {
        row_env <- est_idx[i]; p_i <- per[row_env]; r_i <- ram[row_env]
        for (k in 1:K) {
          s_env <- 0.0
          for (m in 1:M) {
            j <- env_ind_mat[m, row_env]
            if (env_mode == "full") {
              s_env <- s_env + kEnvP_val(p_i, m, k, j)
            } else {
              s_env <- s_env + mu_EnvP[m, k, j]
            }
          }
          inner[i, k] <- kbp[k, p_i] + s_env + br[k, r_i] * sigma_r[k]
        }
      }
      out <- matrix(NA_real_, N_est, K)
      sc_est <- scg[est_idx]
      for (i in 1:N_est) out[i, ] <- softplus(inner[i, ] * sc_est[i], 1)
      out
    }
    
    # baseline
    gdp_imp_base <- impute_from_env(data$env_ind)
    G_base <- matrix(NA_real_, N_env, K)
    G_base[est_idx, ] <- gdp_imp_base
    G_base[eco_idx, ] <- gdp_true
    
    # OND-shifted counterfactuals (lag-aware)
    for (b in seq_along(bin_seq)) {
      sbin <- bin_seq[b]
      env_cf <- data$env_ind
      
      if (!is.null(rn)) {
        keep <- rn %in% target_vars
        for (mname in rn[keep]) {
          m <- which(rn == mname)
          idx_rows <- var_masks[[mname]]
          if (length(idx_rows)) {
            new_idx <- pmax(1L, pmin(NP, env_cf[m, idx_rows] + sbin))
            env_cf[m, idx_rows] <- new_idx
          }
        }
      } else {
        # if no rownames, conservatively shift all rows occurring in the impact window
        idx_rows <- est_idx[ month_from_period(per[est_idx]) %in% impact_months ]
        env_cf[, idx_rows] <- pmax(1L, pmin(NP, env_cf[, idx_rows] + sbin))
      }
      
      gdp_imp_cf <- impute_from_env(env_cf)
      G_cf <- G_base
      G_cf[est_idx, ] <- gdp_imp_cf
      
      
      base_sum <- colSums(G_base[impact_mask_env, , drop = FALSE])
      cf_sum   <- colSums(G_cf  [impact_mask_env, , drop = FALSE])
      pct_chg  <- ifelse(base_sum > 0, (cf_sum / base_sum - 1) * 100, 0)
      arr[b, , d]  <- pct_chg
    }
  }
  
  list(pct_draws = arr, bin_seq = bin_seq, sector_names = sec_names)
}

summarize_bins <- function(posterior, prob = 0.90) {
  arr <- posterior$pct_draws
  B <- dim(arr)[1]; K <- dim(arr)[2]; S <- dim(arr)[3]
  out <- vector("list", K)
  for (k in 1:K) {
    M <- matrix(arr[, k, ], nrow = B, ncol = S)
    med <- apply(M, 1, median, na.rm = TRUE)
    lohi <- t(apply(M, 1, function(v) hpdi(v, prob)))
    out[[k]] <- data.frame(
      bin_shift = posterior$bin_seq,
      y_med = med, y_lo = lohi[,1], y_hi = lohi[,2],
      sector = colnames(arr)[k] %||% posterior$sector_names[k],
      row.names = NULL
    )
  }
  do.call(rbind, out)
}

## -------- OVERALL draws (lag-aware; fixed zero-baseline handling) ----------
build_overall_bin_draws <- function(post, data,
                                    bin_seq = -4:4,
                                    target_vars = c("ind_moisture",
                                                    "ind_lag2_rain","ind_lag3_rain",
                                                    "ind_lag4_rain","ind_lag5_rain","ind_lag6_rain"),
                                    draws = 200,
                                    zero_baseline_policy = c("zero","na")) {
  zero_baseline_policy <- match.arg(zero_baseline_policy)
  
  N_env   <- data$N_env
  N_est   <- data$N_est
  N_eco   <- data$N_eco
  K       <- data$K
  M       <- nrow(data$env_ind)
  NP      <- dim(data$DmatEnv)[2]
  est_idx <- data$est_ind
  eco_idx <- data$eco_ind
  per     <- as.integer(data$period)
  ram     <- as.integer(data$ramadan)
  scg     <- as.numeric(data$scale_gdp); scg[!is.finite(scg) | scg <= 0] <- 1
  P1      <- as.numeric(data$P1)
  rn      <- rownames(data$env_ind)
  
  # lag-aware masks and impact window
  lag_info <- make_lag_row_masks(per, est_idx, target_vars, rn)
  var_masks <- lag_info$masks
  impact_months <- lag_info$months_all
  impact_rows_env <- which(month_from_period(per) %in% impact_months)
  impact_mask_env <- seq_len(N_env) %in% impact_rows_env
  
  S_all <- length(post$a)
  use   <- seq_len(min(draws, S_all))
  
  overall <- matrix(NA_real_, nrow = length(bin_seq), ncol = length(use),
                    dimnames = list(bin = as.character(bin_seq), draw = NULL))
  
  # robust extractors
  get_mu_EnvP_draw <- function(post, s, M, K, NP) {
    x <- post$mu_EnvP; dx <- dim(x)
    if (length(dx) == 4L) array(x[s,,, , drop = FALSE], dim = c(M, K, NP))
    else if (length(dx) == 3L) array(x, dim = c(M, K, NP))
    else stop("mu_EnvP dims: ", paste(dx, collapse="x"))
  }
  get_br_draw <- function(post, s, K) {
    x <- post$br; dx <- dim(x)
    if (length(dx) == 3L) array(x[s,, , drop = FALSE], dim = c(K,2))
    else if (length(dx) == 2L && dx[1]==K && dx[2]==2) array(x, dim = c(K,2))
    else stop("br dims: ", paste(dx, collapse="x"))
  }
  get_sigma_r_draw <- function(post, s, K) {
    x <- post$sigma_r; dx <- dim(x)
    if (length(dx) == 2L) as.numeric(x[s, ])
    else if (length(x) == K) as.numeric(x)
    else stop("sigma_r dims unexpected.")
  }
  get_gdp_true_draw <- function(post, s, N_eco, K) {
    x <- post$gdp_true; dx <- dim(x)
    if (length(dx) == 3L) array(x[s,, , drop = FALSE], dim = c(N_eco, K))
    else if (length(dx) == 2L && dx[1]==N_eco && dx[2]==K) array(x, dim = c(N_eco, K))
    else stop("gdp_true dims: ", paste(dx, collapse="x"))
  }
  
  for (d in seq_along(use)) {
    s <- use[d]
    
    kbp <- matrix(NA_real_, K, data$LP)
    for (k in 1:K) {
      Kper <- gp_periodic_cov_R(P1, post$scalesq_b[s,k], post$length_scale_b[s,k], 26.0)
      kbp[k, ] <- drop(chol(Kper) %*% post$kbp_hat[s, k, ])
    }
    
    mu_EnvP  <- get_mu_EnvP_draw(post, s, M, K, NP)
    br       <- get_br_draw(post, s, K)
    sigma_r  <- get_sigma_r_draw(post, s, K)
    gdp_true <- get_gdp_true_draw(post, s, N_eco, K)
    
    impute_from_env <- function(env_ind_mat) {
      inner <- matrix(0, N_est, K)
      for (i in 1:N_est) {
        row_env <- est_idx[i]; p_i <- per[row_env]; r_i <- ram[row_env]
        for (k in 1:K) {
          s_env <- 0.0
          for (m in 1:M) {
            j <- env_ind_mat[m, row_env]
            s_env <- s_env + mu_EnvP[m, k, j]
          }
          inner[i, k] <- kbp[k, p_i] + s_env + br[k, r_i] * sigma_r[k]
        }
      }
      out <- matrix(NA_real_, N_est, K)
      sc_est <- scg[est_idx]
      for (i in 1:N_est) out[i, ] <- softplus(inner[i, ] * sc_est[i], 1)
      out
    }
    
    gdp_imp_base <- impute_from_env(data$env_ind)
    G_base <- matrix(NA_real_, N_env, K)
    G_base[est_idx, ] <- gdp_imp_base
    G_base[eco_idx, ] <- gdp_true
    
    base_total <- sum(rowSums(G_base[impact_mask_env, , drop = FALSE]))
    
    for (b in seq_along(bin_seq)) {
      sbin <- bin_seq[b]
      env_cf <- data$env_ind
      
      if (!is.null(rn)) {
        keep <- rn %in% target_vars
        for (mname in rn[keep]) {
          m <- which(rn == mname)
          idx_rows <- var_masks[[mname]]
          if (length(idx_rows)) {
            env_cf[m, idx_rows] <- pmax(1L, pmin(NP, env_cf[m, idx_rows] + sbin))
          }
        }
      } else {
        idx_rows <- est_idx[ month_from_period(per[est_idx]) %in% impact_months ]
        env_cf[, idx_rows] <- pmax(1L, pmin(NP, env_cf[, idx_rows] + sbin))
      }
      
      gdp_imp_cf <- impute_from_env(env_cf)
      G_cf <- G_base
      G_cf[est_idx, ] <- gdp_imp_cf
      
      cf_total <- sum(rowSums(G_cf[impact_mask_env, , drop = FALSE]))
      
      if (base_total > 0) {
        overall[b, d] <- (cf_total / base_total - 1) * 100
      } else {
        overall[b, d] <- if (zero_baseline_policy == "zero") 0 else NA_real_
      }
    }
  }
  
  list(pct_overall_draws = overall, bin_seq = bin_seq)
}

summarize_overall <- function(overall, prob = 0.90) {
  m  <- apply(overall$pct_overall_draws, 1, mean, na.rm = TRUE)
  lo <- apply(overall$pct_overall_draws, 1, function(z) hpdi(as.numeric(z), prob)[1])
  hi <- apply(overall$pct_overall_draws, 1, function(z) hpdi(as.numeric(z), prob)[2])
  data.frame(bin_shift = overall$bin_seq, mean = m, lo = lo, hi = hi, row.names = NULL)
}

## -------- mm per bin (robust) ----------
compute_mm_per_bin <- function(data, rain_vars = c("ind_lag2_rain","ind_lag3_rain","ind_lag4_rain","ind_lag5_rain","ind_lag6_rain")) {
  rain_idx <- which(rownames(data$env_ind) %in% rain_vars)
  if (!length(rain_idx)) return(1)  # fallback
  deltas <- unlist(lapply(
    rain_idx,
    function(m) data$DmatEnv[m, 2:dim(data$DmatEnv)[2], 1:(dim(data$DmatEnv)[2]-1)]
  ))
  mm <- mean(deltas, na.rm = TRUE)
  if (!is.finite(mm)) 1 else mm
}

## -------- cases mapping ----------
build_cases_summary <- function(summ_overall, baseline_cases = 100, elasticity = -0.6) {
  trans <- function(pct) baseline_cases * (1 + elasticity * pct / 100)
  data.frame(
    bin_shift = summ_overall$bin_shift,
    mean = trans(summ_overall$mean),
    lo   = trans(summ_overall$lo),
    hi   = trans(summ_overall$hi)
  )
}

## -------- BUILD everything ----------
# 1) sector draws & summary
posterior_bins <- build_sector_bin_draws(
  post, data,
  bin_seq     = -4:4,
  target_vars = c("ind_moisture","ind_lag2_rain","ind_lag3_rain",
                  "ind_lag4_rain","ind_lag5_rain","ind_lag6_rain"),
  draws = 200
)
summ <- summarize_bins(posterior_bins, prob = 0.90)

# 2) overall draws & summary (FIXED zero-baseline handling)
overall_bins <- build_overall_bin_draws(
  post, data,
  bin_seq     = posterior_bins$bin_seq,
  target_vars = c("ind_moisture","ind_lag2_rain","ind_lag3_rain",
                  "ind_lag4_rain","ind_lag5_rain","ind_lag6_rain"),
  draws = 200,
  zero_baseline_policy = "zero"
)
summ_overall <- summarize_overall(overall_bins, prob = 0.90)

# 3) mm per +1 bin
mm_per_bin <- compute_mm_per_bin(data)

# 4) cases mapping (set these to your context)
baseline_cases <- 100   # EDIT

# ---- Elasticity posterior from bgdp_mu using fixed baseline ----
# Returns draws and a compact summary.
elasticity_draws <- function(post, data, softplus_alpha = 1, prob = 0.90) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  hpdi <- function(x, prob = 0.90) {
    x <- sort(as.numeric(x)); n <- length(x); m <- max(1, ceiling(prob*n))
    j <- which.min(x[m:n] - x[1:(n-m+1)])
    c(x[j], x[j+m-1])
  }
  softplus <- function(x, alpha = 1) log1p(exp(alpha*x))/alpha
  logistic <- function(x) 1/(1+exp(-x))
  
  # pull draws
  mu_draw <- if ("bgdp_mu" %in% names(post)) post$bgdp_mu else post$mu_gdp
  a_draw  <- post$a
  
  # dims/indices
  N_env   <- data$N_env
  K       <- data$K
  est_ind <- data$est_ind
  eco_ind <- data$eco_ind
  clove_hist <- data$clove_hist %||% NULL
  
  # fixed GDP baseline from posterior means
  gdp_impute_bar <- apply(post$gdp_impute, c(2,3), mean)  # [N_est,K]
  gdp_true_bar   <- apply(post$gdp_true,   c(2,3), mean)  # [N_eco,K]
  
  gm_fixed <- matrix(NA_real_, N_env, K)
  gm_fixed[est_ind, ] <- gdp_impute_bar
  gm_fixed[eco_ind, ] <- gdp_true_bar
  if (!is.null(clove_hist) && K == 10) gm_fixed[est_ind, 10] <- clove_hist
  
  # reference summaries (fixed, not random)
  S_bar   <- mean(rowSums(log1p(gm_fixed)))                  # \overline{S}
  C_bar   <- mean(rowSums(gm_fixed / (1 + gm_fixed)))        # \overline{C}
  a_bar   <- mean(a_draw)
  mu_bar  <- mean(mu_draw)
  
  # link pieces
  eta_ref   <- a_bar + mu_bar * S_bar
  kappa_ref <- C_bar * (logistic(eta_ref) / softplus(eta_ref, softplus_alpha))
  
  # elasticity draws (linear rescale)
  el_draws <- as.numeric(kappa_ref * mu_draw)
  
  # summary
  hi_lo <- hpdi(el_draws, prob = prob)
  out <- list(
    draws   = el_draws,
    summary = c(
      mean = mean(el_draws),
      sd   = sd(el_draws),
      lo   = hi_lo[1],
      hi   = hi_lo[2],
      prob = prob
    ),
    components = list(
      S_bar = S_bar, C_bar = C_bar, a_bar = a_bar, mu_bar = mu_bar,
      eta_ref = eta_ref, kappa_ref = kappa_ref
    )
  )
  class(out) <- c("elasticity_posterior", class(out))
  out
}

elasticity_draws <- elasticity_draws(post, data)

# hpdi helper
hpdi <- function(x, prob = 0.90) {
  x <- sort(as.numeric(x)); n <- length(x); m <- max(1, ceiling(prob*n))
  j <- which.min(x[m:n] - x[1:(n-m+1)])
  c(x[j], x[j+m-1])
}

# ---- build posterior summary for cases from overall %Change draws + elasticity draws ----
build_cases_from_overall_draws <- function(overall_bins, elasticity_draws,
                                           baseline_cases = 1000, prob = 0.90) {
  pct_mat <- overall_bins$pct_overall_draws
  B <- nrow(pct_mat); S <- ncol(pct_mat)
  
  stopifnot(length(elasticity_draws) >= S)
  elast <- as.numeric(elasticity_draws)[seq_len(S)]
  
  cases_draws <- matrix(NA_real_, B, S)
  for (b in 1:B) {
    cases_draws[b, ] <- baseline_cases * (1 + elast * pct_mat[b, ] / 100)
  }
  
  med  <- apply(cases_draws, 1, median, na.rm = TRUE)
  lohi <- t(apply(cases_draws, 1, function(v) hpdi(v, prob)))
  data.frame(
    bin_shift = overall_bins$bin_seq,
    mean = med, lo = lohi[,1], hi = lohi[,2],
    row.names = NULL
  )
}

hpdi <- function(x, prob = 0.9) {
  x <- sort(as.numeric(x)); n <- length(x)
  m <- ceiling(prob * n)
  j <- which.min(x[m:n] - x[1:(n - m + 1)])
  c(x[j], x[j + m - 1])
}

summ_cases <- build_cases_from_overall_draws(
  overall_bins,
  elasticity_draws$draws,
  baseline_cases = 100,   # <-- set this to your baseline
  prob = 0.90
)

build_summary_from_posterior <- function(posterior_bins, mm_per_bin, prob = 0.90) {
  stopifnot(length(dim(posterior_bins$pct_draws)) == 3)
  B <- dim(posterior_bins$pct_draws)[1]
  K <- dim(posterior_bins$pct_draws)[2]
  S <- dim(posterior_bins$pct_draws)[3]
  
  x_mm <- posterior_bins$bin_seq * mm_per_bin
  s_names <- posterior_bins$sector_names
  if (is.null(s_names) || length(s_names) != K)
    s_names <- paste0("sector_", seq_len(K))
  
  out <- vector("list", K)
  for (k in 1:K) {
    M <- matrix(posterior_bins$pct_draws[, k, ], nrow = B, ncol = S)
    med  <- apply(M, 1, median, na.rm = TRUE)
    lohi <- t(apply(M, 1, function(v) hpdi(v, prob)))
    out[[k]] <- data.frame(
      ond   = x_mm,
      y_med = med,
      y_lo  = lohi[, 1],
      y_hi  = lohi[, 2],
      sector = s_names[k],
      row.names = NULL
    )
  }
  do.call(rbind, out)
}

# estimate mm change per bin (as before)
rain_idx <- which(rownames(data$env_ind) %in%
                    c("ind_lag2_rain","ind_lag3_rain",
                      "ind_lag4_rain","ind_lag5_rain","ind_lag6_rain"))
rain_deltas <- unlist(lapply(rain_idx, function(m)
  data$DmatEnv[m, 2:dim(data$DmatEnv)[2], 1:(dim(data$DmatEnv)[2]-1)]
))
mm_per_bin <- mean(rain_deltas, na.rm = TRUE)

# summarize posterior draws into summary_df
summary_df <- build_summary_from_posterior(posterior_bins, mm_per_bin, prob = 0.9)





#############################################################
############### PLOT ########################################



plot_overall_bins <- function(summ_overall, mm_per_bin,
                              ylab = "%Change in OND through June earnings (total)",
                              xlab = "Change in OND rainfall (mm)",
                              prob_label = "90% band") {
  if (!nrow(summ_overall)) { plot.new(); mtext("No aggregate data", 3); return(invisible()) }
  x_real <- summ_overall$bin_shift * mm_per_bin
  yr <- .safe_range(c(summ_overall$lo, summ_overall$hi, 0))
  
  #op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  #par(mar = c(4, 4, 2, 1), mgp = c(2.1, .7, 0), tcl = -0.2, lend = 1, family = "sans")
  
  plot(x_real, summ_overall$mean, type = "n", bty = "n",
       xlab = xlab, ylab = ylab, ylim = c(-13, 13), xaxt = "n", cex.lab = 1.5)
  title("B.", adj = 0)
  axis(1, at = .pretty_at(x_real))
  abline(h = 0, lty = 2, col = col_border)
  if (any(is.finite(summ_overall$lo) & is.finite(summ_overall$hi))) {
    polygon(c(x_real, rev(x_real)),
            c(summ_overall$lo, rev(summ_overall$hi)),
            col = col_light, border = NA)
  }
  if (any(is.finite(summ_overall$mean))) {
    lines(x_real, summ_overall$mean, col = col_teal_d, lwd = 2)
    points(x_real, summ_overall$mean, pch = 16, col = col_teal, cex = 1)
  }
  #mtext(prob_label, side = 3, cex = 0.8, line = -1.0, adj = 0.02)
}

plot_cases_panel <- function(summ_cases, mm_per_bin,
                             ylab = "Predicted number of cases",
                             xlab = "Change in OND rainfall (mm)",
                             label_note = "Mapped from aggregate %Change earnings") {
  if (!nrow(summ_cases)) { plot.new(); mtext("No cases data", 3); return(invisible()) }
  x_real <- summ_cases$bin_shift * mm_per_bin
  yr <- .safe_range(c(summ_cases$lo, summ_cases$hi))
  
  #op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  #par(mar = c(4, 4, 2, 1), mgp = c(2.1, .7, 0), tcl = -0.2, lend = 1, family = "sans")
  
  plot(x_real, summ_cases$mean, type = "n", bty = "n",
       xlab = xlab, ylab = ylab, ylim = c(93, 107), xaxt = "n", cex.lab = 1.5)
  axis(1, at = .pretty_at(x_real))
  title("C.", adj = 0)
  if (any(is.finite(summ_cases$lo) & is.finite(summ_cases$hi))) {
    polygon(c(x_real, rev(x_real)),
            c(summ_cases$lo, rev(summ_cases$hi)),
            col = col_light, border = NA)
  }
  if (any(is.finite(summ_cases$mean))) {
    lines(x_real, summ_cases$mean, col = col_teal_d, lwd = 2)
    points(x_real, summ_cases$mean, pch = 16, col = col_teal, cex = 1)
  }
  abline(h = 100, lty = 2, col = col_border)
  #mtext(label_note, side = 3, cex = 0.8, line = -1.0, adj = 0.02)
}


plot_sector_strips_strict <- function(summary_df,
                                      xlab = "Change OND rainfall (mm)",
                                      ylab = "Predicted % change in earnings") {
  stopifnot(all(c("ond","y_med","y_lo","y_hi","sector") %in% names(summary_df)))
  
  # coerce types (avoid factor/char issues)
  summary_df$ond    <- as.numeric(as.character(summary_df$ond))
  summary_df$y_med  <- as.numeric(as.character(summary_df$y_med))
  summary_df$y_lo   <- as.numeric(as.character(summary_df$y_lo))
  summary_df$y_hi   <- as.numeric(as.character(summary_df$y_hi))
  summary_df$sector <- as.character(summary_df$sector)
  
  sectors <- unique(summary_df$sector_pretty)
  nS <- length(sectors); if (!nS) { plot(0,0,type="n",axes=FALSE,ann=FALSE); text(0,0,"No data"); return(invisible()) }
  #sectors <- (sectors)  # top→bottom order
  
  # shared axes
  xlim <- range(summary_df$ond, na.rm = TRUE); if (!all(is.finite(xlim))) xlim <- c(-1,1)
  ylim <- range(c(summary_df$y_lo, summary_df$y_hi, 0), na.rm = TRUE); if (!all(is.finite(ylim))) ylim <- c(-1,1)
  yr <- diff(ylim); if (!is.finite(yr) || yr == 0) { ylim <- ylim + c(-1,1); yr <- diff(ylim) }
  
  # safe par save/restore (no layout fields)
  keep <- c("mar","mai","mgp","tcl","cex","las","xaxs","yaxs","lend","xpd","family")
  op <- par(keep); on.exit(par(op), add = TRUE)
  
  # normalized strip geometry
  pad_top <- 0.06; pad_bot <- 0.12; gap <- 0.02
  strip_h <- (1 - pad_top - pad_bot - (nS - 1) * gap) / nS
  
  # colors
  line_col   <- "#0D5E59"
  point_col  <- "#1B9088"
  ribbon_col <- "#E8F4F3"
  zero_col   <- "#BCCAC8"
  
  # base plot occupies current panel; no axes
  par(mar = c(3.6, 15, 1.2, 1.0), xaxs = "i", yaxs = "i", las = 1)
  plot(NA, xlim = xlim, ylim = c(0,1), xlab = "", ylab = "", axes = FALSE, bty = "n")
  title("A.", adj = 0)
  # helper: map real y → strip band [y0,y1]
  to_strip <- function(v, y0, y1) y0 + (v - ylim[1]) / yr * (y1 - y0)
  
  # outer y label
  mtext(ylab, side = 3, line = -1.5, las = 1, cex = .8)
  
  for (i in seq_along(sectors)) {
    s <- sectors[i]
    d <- summary_df[summary_df$sector_pretty == s & is.finite(summary_df$ond), , drop = FALSE]
    if (!nrow(d)) next
    o <- order(d$ond); d <- d[o, ]
    
    y0 <- pad_bot + (nS - i) * (strip_h + gap)
    y1 <- y0 + strip_h
    
    # zero line
    abline(h = to_strip(0, y0, y1), lty = 2, col = zero_col)
    
    # ribbon
    ok_rib <- is.finite(d$y_lo) & is.finite(d$y_hi)
    if (any(ok_rib)) {
      polygon(
        x = c(d$ond[ok_rib], rev(d$ond[ok_rib])),
        y = c(to_strip(d$y_lo[ok_rib], y0, y1), rev(to_strip(d$y_hi[ok_rib], y0, y1))),
        col = ribbon_col, border = NA
      )
    }
    
    # median + points
    ok_med <- is.finite(d$y_med)
    if (any(ok_med)) {
      lines(d$ond[ok_med], to_strip(d$y_med[ok_med], y0, y1), col = line_col, lwd = 2.2)
      points(d$ond[ok_med], to_strip(d$y_med[ok_med], y0, y1), pch = 16, col = point_col, cex = 0.9)
    }
    
    # strip label
    mtext(gsub("_"," ", s), side = 2, at = (y0 + y1)/2, line = 2.7, cex = 0.9)
  }
  
  # single bottom x-axis
  xticks <- pretty(xlim, n = 5)
  axis(1, at = xticks, labels = formatC(xticks, format = "f", digits = 0), line = -1.5)
  mtext(xlab, side = 1, line = 1.5, cex = 1)
  
  invisible(NULL)
}

# internal -> pretty
sector_map <- c(
  agriculture              = "Farming",
  animals                  = "Livestock",
  construction             = "Construction",
  forest                   = "Forestry",
  karafu                   = "Cloves",
  manufacturing_and_repair = "Manufacturing",
  marine                   = "Fishing",
  mwani                    = "Seaweed",
  other_services           = "Other Services",
  public_admin             = "Government",
  trade_and_transport      = "Trade and Transport"
)

# apply mapping
summary_df$sector_pretty <- unname(sector_map[as.character(summary_df$sector)])

# sanity: list any unmatched
unmatched <- unique(summary_df$sector[is.na(summary_df$sector_pretty)])
if (length(unmatched)) warning("Unmapped sector codes: ", paste(unmatched, collapse = ", "))

sectors_desired <- c(
  "Other Services","Cloves","Trade and Transport","Government",
  "Seaweed","Fishing","Manufacturing","Forestry","Construction",
  "Livestock","Farming"
)
# reverse for top→bottom (Farming first)
sectors_plot <- rev(sectors_desired)

pdf("figures/predicted.pdf", height = 6, width = 12)

layout(matrix(c(1,2,3), 1))
# A: strips (does not alter layout pointers)
#par(mar = c(5, 20, 1,1))
plot_sector_strips_strict(summary_df)

# B and C: your other panels (these will now render in cells 2 and 3)
par(mar = c(5.2, 4.2, 1,1))
plot_overall_bins(summ_overall, mm_per_bin)
plot_cases_panel(summ_cases, mm_per_bin)
dev.off()



