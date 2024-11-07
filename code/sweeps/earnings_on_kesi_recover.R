#########################################
######### RECOVER WEATHER ON EARNINGS ###
library(readr)
library(abind)
library(rethinking)
library(posterior)
# 
# 
# set_project_wd <- function(folder){
#   user=Sys.info()[[6]]
#   if(user=="jeffrey_andrews") setwd(paste0("C:/Users/jeffrey_andrews/OneDrive/Documents/", folder))
#   else if(user=="Jeff") setwd(paste0("C:/Users/Jeff/OneDrive/Documents/", folder))
#   else if(user == 'jeffr')  setwd(paste0("C:/Users/jeffr/OneDrive/Documents/", folder))
# }
# set_project_wd("Bio_econ")
# # source("code/script/cleaning/getData.R")
# source("~/functions/utility.R")
# # source("code/Bio_econ.R")


parallelized_function<- function(sweep_list = sweep_list, observed_years = 2){

  num_cores <- detectCores()
  
  # Parallel execution using mclapply
    mclapply(1:nrow(sweep_list), function(sweep) {
      price =  sweep_list[sweep,1]
      inspect = sweep_list[sweep,2]
      ending= paste0("_price_",price, "_inspect_", inspect, ".csv")
      
      earnings_full <- as.data.frame(read_csv(paste0("data/sweeps/earnings_on_kesi/abm/earnings_full", ending)))
      weather <- as.data.frame(read_csv(paste0("data/sweeps/earnings_on_kesi/abm/weather", ending)))
      effort <- as.data.frame(read_csv(paste0("data/sweeps/earnings_on_kesi/abm/effort", ending)))
      kesi <- as.data.frame(read_csv(paste0("data/sweeps/earnings_on_kesi/abm/kesi", ending)))
      total_years = observed_years+10
      nperiods = 26
      data=list(N_obs = observed_years*nperiods,
                NP = nperiods,
                M = ncol(weather),
                N_cases = total_years*nperiods)
      
      #############################
      ######## MAKE GP ############
      
      # Now all we need is a function that loops over the weather data constructing the
      # Necessary DMATs and INDS for the model.
      Env_Ind <- NULL
      EnvDmats <- array(NA, dim = c(ncol(weather), 16, 16))
      for(i in 1:ncol(weather)){
        w_cat=makeGPCat(x = weather[,i], y = weather[,i], L = 16, min_buffer = .8, max_buffer = .8)
        EnvDmat <- makeGPDmat(w_cat[[2]])
        Env_Ind<-cbind(Env_Ind, w_cat[[1]])
        EnvDmats[i,,]<-EnvDmat
      }
      
      data$K = ncol(earnings_full)-1 # Note that this is poorly specified because the last column should be droped
      data$y = earnings_full[(nrow(earnings_full)-data$N_obs+1):nrow(earnings_full),1:data$K] *10
      data$N_eco <- nrow(data$y)
      data$N_env <- data$N
      data$N_est = data$N_cases- data$N_obs
      data$N_overlap = data$N_obs
      data$sver <- sample(1:2, data$N_obs, replace = T)
      data$NV <- length(unique(data$sver))
      data$N = total_years*nperiods
      data$LP <- nperiods
      data$env_ind = t(Env_Ind[(nrow(earnings_full)-data$N+1):nrow(earnings_full),])
      data$N_env <- data$N
      # Inds
      data$est_ind <- 1:data$N_est
      data$obs_ind <- (data$N_est+1):data$N
      data$eco_ind <- data$obs_ind
      # other
      data$DmatEnv <- EnvDmats/10
      data$scale_gdp <- rep(1, data$N_obs+data$N_est)
      data$P1 <- 1:nperiods
      softplus<-1
      data$L = total_years
      data$period <- rep(data$P1, data$L)
      data$NP <- ncol(data$DmatEnv[,,1])
      data$ramadan<-sample(1:2, data$N, replace = T)
      data$softplus_alpha = 1
      period_means <- matrix(NA, ncol = data$K, nrow = data$LP)
      for(i in 1:data$LP){
        for(k in 1:data$K){
          inds = which(data$period[data$obs_ind] == i)
          means = data$y[inds,k]
          weights = rep(1/length(inds), length(inds))
          sum = sum(weights)
          # period_means[i, k] = sum(means*weights)/sum
           period_means[i, k] <- mean(data$y[which(data$period[data$obs_ind] == i),k]) # Note this this the unweighted Version
        }
      }
      
      data$period_means = period_means-colMeans(period_means)
      # data$period_means = period_means
      data$k_means = colMeans(period_means)
      data$y2 <-earnings_full[(nrow(earnings_full)-data$N_env+1):nrow(earnings_full),1:data$K] *10
      data$kesi <- kesi[(nrow(kesi)-data$N_env+1):nrow(kesi),1]
      data$denv_ind <- matrix(NA, nrow = 1, ncol = length(data$env_ind[1,]))
      data$denv_ind[1,] <- data$env_ind[1,]
      data$DmatDEnv <- array(NA, c(1, data$NP, data$NP ))
      data$DmatDEnv[1,,] <- data$DmatEnv[1,,] 
      data$DM =1 
      data$scale_gdp <- rep(1, data$N_obs+data$N_est)
      data$timber_prices <- rnorm(data$N_obs+data$N_est)
      # Make matrix for DMAT
      L = total_years
      Dmat = matrix(NA, L, L)
      for(i in 1:L){
        for(j in 1:L){
          Dmat[i, j] = abs(i-j)
        }  
      }
      data$DmatX= Dmat
      data$year = rep(1:total_years, each= nperiods)
      data$target_var = .05
      
      # First Check to see if we can see the full set of earnings what the parameter values are
      
      data$y <- data$y +.01
      data$y2 <- data$y2 +.01
      
      ###########################
      ##### FULL DATA ###########
      # Second see what happens if the only see the last segrement on of the data!
      mod<-cmdstanr::cmdstan_model("code/stan_models/recover_earnings_on_kesi_without_imputation.stan")
      #mod2<-cmdstanr::cmdstan_model("C:/Users/jeffr/OneDrive/Documents/Bio_econ/code/stan_models/test_on_multivariate_normal.stan")
      pf <- mod$pathfinder(data = data)
      #pf2 <- mod2$pathfinder(data = data)
      init_list<-get_init_list(pf)
      a = mod$sample(parallel_chains =4,
                     chains =4,
                     data = data,
                     init = list(init_list, init_list, init_list, init_list),
                     iter_warmup =500, iter_sampling = 500, refresh = 50)
      post<-extract.samples2(a)
      output1_full<- (post$bgdp[,1]*post$sigma_gdp+post$mu_gdp)
      output2_full<- (post$bgdp[,2]*post$sigma_gdp+post$mu_gdp)
      #######################
      ##### IMPUTED DATA ####
      
      # Second see what happens if the only see the last segrement on of the data!
      mod<-cmdstanr::cmdstan_model("code/stan_models/recover_earnings_on_kesi.stan")
      # mod2<-cmdstanr::cmdstan_model("C:/Users/jeffr/OneDrive/Documents/Bio_econ/code/stan_models/test_on_multivariate_normal.stan")
      
      tryCatch({
        pf <- mod$pathfinder(data = data)
      }, error = function(e) {
        # Log the error details
        message("An error occurred in sweep,", ending, ". The pathfinder function returns the following error: ", e$message,)
        # Return a default value or handle the error as needed
      })

      # pf2 <- mod2$pathfinder(data = data)
      init_list<-get_init_list(pf)
      
      a2 = mod$sample(parallel_chains =4,
                     chains =4,
                     data = data,
                     init = list(init_list, init_list, init_list, init_list),
                     iter_warmup =500, iter_sampling = 500, refresh = 10)
      
      post2<-extract.samples2(a2)
      output1_imputed<- (post2$bgdp[,1]*post2$sigma_gdp+post2$mu_gdp)
      output2_imputed<- (post2$bgdp[,2]*post2$sigma_gdp+post2$mu_gdp)
  
      
      write.csv(output1_full, paste0("data/sweeps/earnings_on_kesi/stan/good1_observed_years", observed_years, ending))
      write.csv(output1_imputed, paste0("data/sweeps/earnings_on_kesi/stan/good1_imputed_years", observed_years, ending))
      
      write.csv(output2_full, paste0("data/sweeps/earnings_on_kesi/stan/good2_observed_years", observed_years, ending))
      write.csv(output2_imputed, paste0("data/sweeps/earnings_on_kesi/stan/good2_imputed_years", observed_years, ending))
          
      }, mc.cores = num_cores)
}# End Function


sweep_list <- as.data.frame(read_csv("data/sweeps/earnings_on_kesi/sweep_list.csv"))

parallelized_function(sweep_list = sweep_list)

