

# Main function
parallelized_function <- function(sweep_list, observed_years = 10) {
  
  # Set up the number of cores to use
  num_cores <- detectCores()
  
  # Parallel execution using mclapply
  mclapply(1:nrow(sweep_list), function(sweep) {
    
    
    
    gamma11 <- sweep_list[sweep, 1]
    wage_noise <- sweep_list[sweep, 2]
    ending <- paste0("_noise_", wage_noise, "_gamma_", gamma11, ".csv")
    
    earnings <- as.data.frame(read_csv(paste0("data/sweeps/weather_on_earnings/abm/earnings_observed", ending)))
    earnings_full <- as.data.frame(read_csv(paste0("data/sweeps/weather_on_earnings/abm/earnings_full_observed", ending)))
    earnings_counterfactual <- as.data.frame(read_csv(paste0("data/sweeps/weather_on_earnings/abm/earnings_counterfactual", ending)))
    earnings_full_counterfactual <- as.data.frame(read_csv(paste0("data/sweeps/weather_on_earnings/abm/earnings_full_counterfactual", ending)))
    weather <- as.data.frame(read_csv(paste0("data/sweeps/weather_on_earnings/abm/weather_observed", ending)))
    weather_counterfactual <- as.data.frame(read_csv(paste0("data/sweeps/weather_on_earnings/abm/weather_counterfactual", ending)))
    
    total_years <- observed_years + 1
    nperiods <- 26
    data <- list(N_obs = observed_years * nperiods,
                 NP = nperiods,
                 M = ncol(weather),
                 N_cases = total_years * nperiods)
    
    data$K <- ncol(earnings_full)-1
    data$N_est <- data$N_cases - data$N_obs
    data$N <- total_years * nperiods
    data$y <- earnings_full[(nrow(earnings_full) - data$N_obs + 1):nrow(earnings_full), 1:data$K] * 10
    data$y[data$y<0] <- 0
    
    data$N_eco <- nrow(data$y)
    data$N_env <- data$N
    
    
    Env_Ind <- NULL
    EnvDmats <- array(NA, dim = c(ncol(weather), 16, 16))
    for(j in 1:ncol(weather)) {
      w_cat <- makeGPCat(x = weather[(nrow(weather) - data$N_env + 1):nrow(weather), j], y = weather_counterfactual[(nrow(weather) - data$N_env + 1):nrow(weather), j],
                         L = 16, min_buffer = .8, max_buffer =  1)
      EnvDmat <- makeGPDmat(w_cat[[2]])
      Env_Ind <- cbind(Env_Ind, w_cat[[1]])
      EnvDmats[j, , ] <- EnvDmat
    }
    
    (nrow(weather) - data$N_obs + 1):nrow(weather)
    
    Env_Ind_Counter <- NULL
    for(j in 1:ncol(weather)) {
      w_cat <- makeGPCat(x = weather_counterfactual[(nrow(weather) - data$N_env + 1):nrow(weather), j], y = weather[(nrow(weather) - data$N_env + 1):nrow(weather), j],
                         L = 16, min_buffer = .8, max_buffer =  1)
      Env_Ind_Counter <- cbind(Env_Ind_Counter, w_cat[[1]])
    }
    
    data$N_overlap <- data$N_obs
    data$sver <- sample(1:2, data$N_obs, replace = TRUE)
    data$NV <- length(unique(data$sver))
    data$LP <- nperiods
    data$env_ind <- t(Env_Ind[, ])
    data$env_ind_counter <- t(Env_Ind_Counter[, ])
    data$y_counterfactual <- earnings_full_counterfactual[(nrow(earnings_full_counterfactual) - data$N_obs + 1):nrow(earnings_full_counterfactual), 1:data$K] * 10
    data$y_counterfactual[data$y_counterfactual<0] <- 0
    data$N_env <- data$N
    data$est_ind <- 1:data$N_est
    data$obs_ind <- (data$N_est + 1):data$N
    data$eco_ind <- data$obs_ind
    data$DmatEnv <- EnvDmats / 10
    data$scale_gdp <- rep(1, data$N_obs + data$N_est)
    data$P1 <- 1:nperiods
    softplus <- 1
    data$L <- total_years
    data$period <- rep(data$P1, data$L)
    data$NP <- ncol(data$DmatEnv[, , 1])
    data$ramadan <- sample(1:2, data$N, replace = TRUE)
    data$softplus_alpha <- 1
    data$target_var = .01
    
    period_means <- matrix(NA, ncol = data$K, nrow = data$LP)
    for(p in 1:data$LP) {
      for(k in 1:data$K) {
        inds <- which(data$period[data$obs_ind] == p)
        means <- data$y[inds, k]
        weights <- rep(1 / length(inds), length(inds))
        sum_weights <- sum(weights)
        period_means[p, k] <- sum(means * weights) / sum_weights
      }
    }
    
    data$period_means = period_means-colMeans(period_means)
    #data$period_means <- period_means 
    data$k_means <- colMeans(period_means)
    data$y = data$y +.01

    mod <- cmdstanr::cmdstan_model("code/stan_models/recover_env_on_earnings.stan")
    pf <- mod$pathfinder(data = data)
    init_list <- get_init_list(pf)
    
    a <- mod$sample(parallel_chains = 4,
                    chains = 4,
                    data = data,
                    init = list(init_list, init_list, init_list, init_list),
                    iter_warmup = 500, iter_sampling = 500, refresh = 10)
    
    post <- extract.samples2(a)
    # 
    # nsims <- nrow(post$br[, , 1])
    # for(k in 1:data$K) {
    #   plot(rep(1, nsims) + rnorm(nsims, 0, .1), post$gdp_true[, 1, k],  xlim = c(1, 27), 
    #        col = col.alpha("#249EA0", .1), ylim = c(0, 30), pch = 16,
    #        xlab = "Period", ylab = "Earnings")
    #   for(p in 2:26) {
    #     points(rep(p, nsims) + rnorm(nsims, 0, .1), post$gdp_true[, p, k],  xlim = c(1, 27),
    #            col = col.alpha("#249EA0", .1), ylim = c(0, 30),
    #            xlab = "Period", ylab = "Earnings")
    #   }
    #   lines(1:26, data$y[1:26, k], type = "l", ylim = c(0, 30), lwd = 3, col = "#FAAB36")
    # }
    # 
    # No idea why it cannot ever find this function
    softplus <- function(x, alpha) log1p(exp(x*alpha))/alpha
    
    
    nsims = 1200
    mu <- array(NA, dim = c(nsims, data$N_eco, data$K))
    for(e in 1:data$N_eco) {
      for(k in 1:data$K) {
        environment_effectP <- matrix(rep(NA, data$M * nrow(post$br[, , 1])), ncol = data$M)
        for(m in 1:data$M) {
          environment_effectP[, m] <- with(post, kEnvP[, data$period[data$eco_ind[e]], m, k, data$env_ind[m, data$eco_ind[e]]])
        }
        mu[, e, k] <- with(post, softplus(gamma_mu[, k] +
                                          + kbp[, k, data$period[data$eco_ind[e]]] 
                                          + rowSums(environment_effectP)
                                          + br[, k, data$ramadan[data$eco_ind[e]]] * sigma_r[, k] + b_ver[, k, data$sver[e]],
                                          data$softplus_alpha))
      }
    }
    
    output1 <- array(NA, dim = c(nsims, data$N_eco, data$K))
    for(e in 1:data$N_eco) {
      for(k in 1:data$K) {
        alpha <- mu[, e, k] * post$phi[, k]
        output1[, e, k] <- rgamma(nsims, alpha, post$phi[, k])
      }
    }
    
    mu <- array(NA, dim = c(nsims, data$N_eco, data$K))
    for(e in 1:data$N_eco) {
      for(k in 1:data$K) {
        environment_effectP <- matrix(rep(NA, data$M * nrow(post$br[, , 1])), ncol = data$M)
        for(m in 1:data$M) {
          environment_effectP[, m] <- with(post, kEnvP[, data$period[data$eco_ind[e]], m, k, data$env_ind_counter[m, data$eco_ind[e]]])
        }
        mu[, e, k] <- with(post, softplus(gamma_mu[, k] +
                                          + kbp[, k, data$period[data$eco_ind[e]]] 
                                          + rowSums(environment_effectP)
                                          + br[, k, data$ramadan[data$eco_ind[e]]] * sigma_r[, k] + b_ver[, k, data$sver[e]],
                                          data$softplus_alpha))
      }
    }
    
    output2 <- array(NA, dim = c(nsims, data$N_eco, data$K))
    for(e in 1:data$N_eco) {
      for(k in 1:data$K) {
        alpha <- mu[, e, k] * post$phi[, k]
        output2[, e, k] <- rgamma(nsims, alpha, post$phi[, k])
      }
    }
    
    write.csv(data$y_counterfactual[, 1], paste0("data/sweeps/weather_on_earnings/stan/y_counterfactual_1_observed_years_", observed_years, ending))
    write.csv(data$y_counterfactual[, 2], paste0("data/sweeps/weather_on_earnings/stan/y_counterfactual_2_observed_years_", observed_years, ending))
    write.csv(output2[, ,1], paste0("data/sweeps/weather_on_earnings/stan/y_estimated_1_observed_years_", observed_years, ending))
    write.csv(output2[, ,2], paste0("data/sweeps/weather_on_earnings/stan/y_estimated_2_observed_years_", observed_years, ending))
    
  }, mc.cores = num_cores)
  
} # End function


sweep_list <- as.data.frame(read_csv("data/sweeps/weather_on_earnings/sweep_list.csv"))


parallelized_function(sweep_list = sweep_list)

