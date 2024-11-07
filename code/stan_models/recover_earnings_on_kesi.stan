functions {
  // Smoothed Exponential Gaussian Process
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];  // Get the dimension of the input matrix x
    matrix[N, N] K;  // Initialize the covariance matrix K
    for (i in 1 : (N - 1)) {  // Loop over rows of the matrix, except the last
      K[i, i] = sq_alpha + delta;  // Set the diagonal elements
      for (j in (i + 1) : N) {  // Loop over columns of the matrix
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i, j]));  // Compute the off-diagonal elements using the Gaussian kernel
        K[j, i] = K[i, j];  // Ensure the matrix is symmetric
      }
    }
    K[N, N] = sq_alpha + delta;  // Set the last diagonal element
    return K;  // Return the covariance matrix
  }
  
  // Period Gaussian Process
  matrix cov_periodic(matrix x, real sq_alpha, real sq_rho, real delta, real period) {
    int NOBS = dims(x)[1];  // Get the number of observations
    matrix[NOBS, NOBS] K;  // Initialize the covariance matrix K
    for (i in 1 : (NOBS - 1)) {  // Loop over rows of the matrix, except the last
      K[i, i] = sq_alpha + delta;  // Set the diagonal elements
      for (j in (i + 1) : NOBS) {  // Loop over columns of the matrix
        K[i, j] = (sq_alpha ^ 2) * exp(-2 * (sin(pi() * x[i, j] / period) ^ 2) / (sq_rho ^ 2));  // Compute the off-diagonal elements using the periodic kernel
        K[j, i] = K[i, j];  // Ensure the matrix is symmetric
      }
    }
    K[NOBS, NOBS] = sq_alpha + delta;  // Set the last diagonal element
    return K;  // Return the covariance matrix
  }
  
  // Merge missing for imputation 
  real softplus(real x, real alpha) {
    real output;  // Declare the output variable
    output = log1p(exp(x * alpha)) / alpha;  // Compute the softplus function
    return output;  // Return the result
  }
}

data {
  // For GDP Estimation
  int<lower=0> N_cases;  // Number of Observed Cases for Prediction
  int<lower=0> N_obs;  // Number of time points where cases are observed and we have comparable income data
  int<lower=0> N_est;  // Number of time points GDP data is unavailable but cases are
  int<lower=0> N_eco;  // Number of economic observations
  int<lower=0> N_env;  // Number of environmental observations
  int<lower=0> K;  // Number of sectors
  int<lower=0> NP;  // Number of Percentiles for GPs
  int M;  // Number of Environmental Predictors
  int DM;  // Number of Environmental Predictors having a direct effect on Kesi
  int NV;  // Number of Survey Versions
  int LP;  // Number of Periods
  array[N_obs] int obs_ind;  // A list of indexes for observed data
  array[N_est] int est_ind;  // A list of indexes for data where GDP is to be estimated
  array[N_eco] int eco_ind;  // A list of indexes for economic data
  array[M, N_env] int env_ind;  // Indexed value for environmental predictors of economy
  array[M] matrix[NP, NP] DmatEnv;  // Array of Distance Matrices for Environmental Predictors
  vector[N_env] scale_gdp;  // Scaling factors for GDP
  array[N_eco] int sver;  // GDP per sector
  array[N_env] int period;  // Indicator for particular month
  matrix[N_eco, K] y;  // GDP per sector
  
  // For Kesi Predictions
  int<lower=0> L;  // Number of Years 
  array[N_env] int ramadan;  // Ramadan indicator
  array[LP] real P1;  // Total number of periods
  matrix[LP, K] period_means;  // Mean values for each period and sector
  vector[K] k_means;  // Mean values for each sector
  real softplus_alpha;  // Parameter for the softplus function
  array[DM, N_env] int denv_ind;  // Indexed value for environmental predictors of Kesi
  array[DM] matrix[NP, NP] DmatDEnv;  // Array of Distance Matrices for Direct Environmental Predictors
  vector[N_env] timber_prices;  // Prices of timber
  matrix[L, L] DmatX;  // Year distance matrix
  array[N_env] int year;  // Time period ID 1:N_est
  array[N_cases] int kesi;  // Number of observed kesi in the time period
  real target_var;  // Target variance
}

parameters {
  matrix<lower=0>[N_est, K] gdp_impute;  // Imputed GDP values
  matrix<lower=0>[N_eco, K] gdp_true;  // True GDP values
  array[K] vector[LP] kbp_hat;  // Raw periodic effects for each sector

  // Period Effect
  array[K] real<lower=0> length_scale_b;  // Length scale for periodic GP in Gamma likelihood
  array[K] real<lower=0> scalesq_b;  // Scale parameter for periodic GP in Gamma likelihood
  
  // PERIOD SPECIFIC ENVIRONMENTAL EFFECTS ON KESI 
  array[M, K] real<lower=0> etasq_EnvP;  // Environmental effect scale parameters
  array[M, K] real<lower=0> rhosq_EnvP;  // Environmental effect length scales
  array[M, K] real<lower=0> sigma_EnvP;  // Environmental effect noise terms
  array[LP, M, K] vector[NP] zEnvP;  // Latent variables for environmental effects
  array[M, K] real<lower=0> cov_EnvP;  // Covariance parameters for environmental effects
  array[M, K] vector[NP] mu_EnvP;  // Mean parameters for environmental effects
  
  array[K] real<lower=0> sigma_r;  // Standard deviation for Ramadan effect
  array[K] vector[2] br;  // Ramadan effect parameters
  array[K] real<lower=0> phi;  // Parameter for Gamma distribution
  
  matrix[NV, K] z_ver;  // Raw random effect for household survey version 
  cholesky_factor_corr[NV] L_ver;  // Cholesky factor for Survey Version
  vector<lower=0>[NV] Sigma_ver;  // Standard deviations for survey versions
  vector[K] gamma_mu;  // Mean parameters for gamma distribution
  
  // Kesi Prediction
  real a;  // Intercept
  real btimber;  // Timber price effect
  vector[K] bgdp;  // GDP effect parameters
  real<lower=0> sigma_r_k;  // Standard deviation for Kesi Ramadan effect
  vector[2] br_k;  // Ramadan effect parameters for Kesi
  real<lower=0> sigmaX;  // Variance parameter for periodic GP in Poisson likelihood
  real<lower=0> etasq;  // Scale parameter for periodic GP in Poisson likelihood
  vector[NV] b_sver_kesi;  // Survey version effect for Kesi
  real<lower = 0> sigma_sver_kesi;  // Standard deviation for survey version effect
  real<lower=0> omega;  // Overdispersion parameter for Kesi
  real<lower=0> sigma_gdp;  // Standard deviation for GDP
  real mu_gdp;  // Mean parameter for GDP
  array[DM] real<lower=0> etasq_DEnv;  // Scale parameters for direct environmental effects
  array[DM] real<lower=0> rhosq_DEnv;  // Length scales for direct environmental effects
  array[DM] real<lower=0> sigma_DEnv;  // Noise terms for direct environmental effects
  array[DM] vector[NP] zDEnv;  // Latent variables for direct environmental effects
  
  matrix[K, N_eco] zP;  // Latent variables for periodic GP
  real<lower=0> rhosq;  // Covariance parameter for periodic GP in Poisson likelihood
  real<lower=0> length_scale;  // Length scale for periodic GP in Poisson likelihood
  real<lower=0> scalesq;  // Scale parameter for periodic GP in Poisson likelihood
  vector[L] zX;  // Latent variables for year effect
  vector[LP] kp_hat;  // Raw periodic effects for Kesi
  

}

transformed parameters {
  matrix[N_est, K] gdp_impute_mu;  // Mean imputed GDP values
  matrix[N_eco, K] gdp_error;  // Error in GDP values
  matrix[K, NV] b_ver;  // Random effects for survey versions
  
  // Periodic Effect of Sectors
  array[K] vector[LP] kbp;  // Periodic effects for each sector
  array[K] cov_matrix[LP] SIGMABP;  // Covariance matrices for periodic effects
  
  // Period Specific Effect of ENV on Sectors
  array[LP, M, K] vector[NP] kEnvP;  // Environmental effects for each period and sector
  array[M, K] matrix[NP, NP] LEnv_sigmaP;  // Cholesky factors for environmental effects
  array[M, K] matrix[NP, NP] SIGMAEnvP;  // Covariance matrices for environmental effects
  
  // Direct Effect on Kesi 
  array[DM] vector[NP] kDEnv;  // Direct environmental effects
  array[DM] matrix[NP, NP] LDEnv_sigma;  // Cholesky factors for direct environmental effects
  array[DM] matrix[NP, NP] SIGMADEnv;  // Covariance matrices for direct environmental effects
  
  cov_matrix[LP] SIGMAP;  // Covariance matrix for periodic GP
  vector[L] kX;  // Year effect
  matrix[L, L] LX_sigma;  // Cholesky factor for year effect
  matrix[L, L] SIGMAX;  // Covariance matrix for year effect
  vector[LP] kp;  // Periodic effects for Kesi

  // Environment on Sector Income by period
  for (m in 1:M) {
    for (k in 1 : K) {
      for (p in 1:LP) {
        SIGMAEnvP[m, k] = cov_GPL2(DmatEnv[m], etasq_EnvP[m, k], rhosq_EnvP[m, k], sigma_EnvP[m, k]);  // Compute the covariance matrix
        LEnv_sigmaP[m, k] = cholesky_decompose(SIGMAEnvP[m, k]);  // Compute the Cholesky factor
        kEnvP[p, m, k] = LEnv_sigmaP[m, k] * (zEnvP[p, m, k] * cov_EnvP[m, k] + mu_EnvP[m, k]);  // Compute the environmental effects
      }
    }
  }
  
  b_ver = (diag_pre_multiply(Sigma_ver, L_ver) * z_ver)';  // Compute the random effects for survey versions
  
  // Periodic effects for Sectors
  for (i in 1 : K) {
    SIGMABP[i] = gp_periodic_cov(P1, scalesq_b[i], length_scale_b[i], 26.0);  // Compute the covariance matrix
    kbp[i] = cholesky_decompose(SIGMABP[i]) * kbp_hat[i];  // Compute the periodic effects
  }
  
  // Impute GDP values based on environmental predictors and training data parameters
  for (i in 1 : N_est) {
    for (k in 1 : K) {
      vector[M] environment_effectP;
      for (m in 1:M) {
        environment_effectP[m] = kEnvP[period[est_ind[i]], m, k][env_ind[m, est_ind[i]]];  // Compute the environmental effects
      }
      gdp_impute_mu[i, k] = softplus(
        (kbp[k][period[est_ind[i]]]
        + sum(environment_effectP)
        + br[k][ramadan[est_ind[i]]] * sigma_r[k])
        * scale_gdp[est_ind[i]], softplus_alpha);  // Compute the imputed GDP values
    }
  }
  
  // Compute the errors in GDP values due to survey versions
  for (i in 1 : N_eco) {
    for (k in 1 : K) {
      gdp_error[i, k] = softplus(y[i, k] - b_ver[k, sver[i]], softplus_alpha);  // Compute the error in GDP values
    }
  }
  
  // Periodic effects for Kesi
  SIGMAP = gp_periodic_cov(P1, scalesq, length_scale, 26.0);  // Compute the covariance matrix
  kp = cholesky_decompose(SIGMAP) * kp_hat;  // Compute the periodic effects
  
  // Compose GPs for year effect
  SIGMAX = cov_GPL2(DmatX, etasq, rhosq, sigmaX);  // Compute the covariance matrix
  LX_sigma = cholesky_decompose(SIGMAX);  // Compute the Cholesky factor
  kX = LX_sigma * zX;  // Compute the year effect
  
  // Direct effect of Environment on Cases
  for (m in 1:DM) {
    SIGMADEnv[m] = cov_GPL2(DmatDEnv[m], etasq_DEnv[m], rhosq_DEnv[m], sigma_DEnv[m]);  // Compute the covariance matrix
    LDEnv_sigma[m] = cholesky_decompose(SIGMADEnv[m]);  // Compute the Cholesky factor
    kDEnv[m] = LDEnv_sigma[m] * zDEnv[m];  // Compute the direct environmental effects
  }
}

model {
  matrix[N_eco, K] mu;  // Mean GDP values
  matrix[N_env, K] gdp_merged;  // Merged GDP values
  vector[N_cases] lambda;  // Lambda parameter for Kesi

  // Priors for periodic GP parameters in Gamma likelihood
  to_vector(length_scale_b) ~ normal(0, .2);
  to_vector(scalesq_b) ~ normal(1, 1);
  
  // Priors for gamma_mu
  for (k in 1:K) gamma_mu[k] ~ normal(k_means[k], 1);
  
  // Priors for periodic effects
  for (i in 1:LP) {
    for (k in 1:K) {
      kbp_hat[k][i] ~  normal(period_means[i,k], .1);
    }
  }
  
  // Priors for Ramadan effects
  for (i in 1:K) {
    br[i] ~ normal(0, 3);
    sigma_r[i] ~ exponential(1);
  }
  
  // Priors for Gamma scale parameter
  for (i in 1:K) phi[i] ~ normal(mean(y[, i]) / variance(y[, i]) + 1.5, .3);
  
  // Priors for survey version effects
  L_ver ~ lkj_corr_cholesky(5);
  Sigma_ver ~ exponential(1);
  to_vector(z_ver) ~ normal(0, 1);
  
  // Priors for environmental effects on Kesi
  for (m in 1:M) {
    for (k in 1:K) {
      rhosq_EnvP[m, k] ~ exponential(1);
      etasq_EnvP[m, k] ~ exponential(1);
      sigma_EnvP[m, k] ~ exponential(1);
      cov_EnvP[m, k] ~ exponential(1);
      mu_EnvP[m, k] ~ normal(0, 1);
      for (p in 1:LP) {
        to_vector(zEnvP[p, m, k]) ~ normal(0, 1);
      }
    }
  }
  
  // Priors for Kesi
  to_vector(kp_hat) ~ normal(0, 1);
  to_vector(zP) ~ normal(0, 1);
  rhosq ~ exponential(1);
  etasq ~ exponential(1);
  sigmaX ~ exponential(1);
  
  // Priors for periodic GP parameters
  length_scale ~ normal(0, .1);
  scalesq ~ normal(1, 1);
  to_vector(zX) ~ normal(0, 1);
  
  // Priors for direct environmental effects
  for (m in 1:DM) {
    rhosq_DEnv[m] ~ exponential(1);
    etasq_DEnv[m] ~ exponential(1);
    sigma_DEnv[m] ~ exponential(1);
    to_vector(zDEnv[m]) ~ normal(0, 1);
  }
  
  // Priors for Kesi predictors
  br_k ~ normal(0, 1);
  sigma_r_k ~ exponential(1);
  btimber ~ normal(0, 1);
  b_sver_kesi ~ normal(0, 1);
  sigma_sver_kesi ~ exponential(1);
  a ~ normal(mean(kesi), 1 );
  omega ~ exponential(1);
  sigma_gdp ~ exponential(1);
  mu_gdp ~ normal(0, 5);
  for (i in 1:K) bgdp[i] ~ normal(0, 5);
  
  // Estimate parameters for the impact of environmental variables on each economic sector
  for (i in 1:N_eco) {
    for (k in 1:K) {
      vector[M] environment_effectP;
      for (m in 1:M) {
        environment_effectP[m] = kEnvP[period[eco_ind[i]], m, k][env_ind[m, eco_ind[i]]];
      }
      mu[i, k] = softplus(
                       kbp[k][period[eco_ind[i]]] 
                     + sum(environment_effectP)
                     + br[k][ramadan[eco_ind[i]]] * sigma_r[k] + b_ver[k, sver[i]],
                     softplus_alpha);  // Compute the mean GDP values
                     //print(mu[i, k]);
    }
  }
  
  // Likelihood for observed GDP values
  for (i in 1:K) {
    vector[N_eco] alpha;
    alpha = mu[:, i] * phi[i];
    //print(alpha);
    y[:, i] ~ gamma(alpha, phi[i]);
  }
  
  // Likelihood for imputed GDP values
  for (i in 1:K) {
    gdp_impute[:, i] ~ gamma(gdp_impute_mu[:, i] * phi[i], phi[i]);
  }
  // 
  // Deal with measurement error caused by survey versions
  for (t in 1:N_eco) {
    for (k in 1:K) {
      real wanted_var = sqrt((y[t, k] * target_var) / y[t, k]);
      gdp_true[t, k] ~ gamma(gdp_error[t, k] * wanted_var, wanted_var);
    }
  }

  // Assign imputed values to merged GDP
  for (i in 1:N_est) {
    for (k in 1:K) {
      gdp_merged[i, k] = gdp_impute[i, k];
    }
  }

  // Assign real values to merged GDP
  for (i in 1:N_eco) {
    for (k in 1:K) {
      gdp_merged[eco_ind[i], k] = gdp_true[i, k];
    }
  }

  // Compute lambda for Kesi
  for (i in 1:N_cases) {
    vector[DM] direct_environmental_effect;
    for (m in 1:DM) {
      direct_environmental_effect[m] = kDEnv[m][denv_ind[m, i]];
    }

    lambda[i] = softplus(a +
      kp[period[i]] +
      dot_product((bgdp * sigma_gdp + mu_gdp), log1p(to_vector(gdp_merged[i,]))) +
      kX[year[i]] +
      br_k[ramadan[i]] * sigma_r_k +
      sum(direct_environmental_effect) +
      btimber * timber_prices[i], softplus_alpha);  // Compute lambda for Kesi
  }

  // Likelihood for observed Kesi
  for (i in 1:N_cases) {
    kesi[i] ~ neg_binomial_2(lambda[i], omega);
  }
}


// functions {
//   // Smoothed Exponential gaussian Process
//   matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
//     int N = dims(x)[1];
//     matrix[N, N] K;
//     for (i in 1 : (N - 1)) {
//       K[i, i] = sq_alpha + delta;
//       for (j in (i + 1) : N) {
//         K[i, j] = sq_alpha * exp(-sq_rho * square(x[i, j]));
//         K[j, i] = K[i, j];
//       }
//     }
//     K[N, N] = sq_alpha + delta;
//     return K;
//   }
//   // Period Gaussian Process
//   matrix cov_periodic(matrix x, real sq_alpha, real sq_rho, real delta,
//                       real period) {
//     int NOBS = dims(x)[1];
//     matrix[NOBS, NOBS] K;
//     for (i in 1 : (NOBS - 1)) {
//       K[i, i] = sq_alpha + delta;
//       for (j in (i + 1) : NOBS) {
//         K[i, j] = (sq_alpha ^ 2)
//                   * exp(-2 * (sin(pi() * x[i, j] / period) ^ 2)
//                         / (sq_rho ^ 2));
//         K[j, i] = K[i, j];
//       }
//     }
//     K[NOBS, NOBS] = sq_alpha + delta;
//     return K;
//   }
//   // Merge missing for imputation 
//   real softplus(real x, real alpha) {
//     real output;
//     output =log1p(exp(x*alpha))/alpha; 
//     return output;
//   }
// }
// data {
//   // For GDP Estimation
//   int<lower=0> N_cases; // Number of Observed Cases for Prediction
//   int<lower=0> N_obs; // Number of time where cases are obsereved and we have comparable income data
//   int<lower=0> N_est; // Number of time points gdp data is unavailible but cases are.
//   int<lower=0> N_eco; // Number of economic obervations
//   int<lower=0> N_env; // Number of economic obervations
//   int<lower=0> K; // Number of sectors
//   int<lower=0> NP; // Number of Percentiles for GPs. 
//   int M; // Number of Environmental Predictors
//   int DM; // Number of Environmental Predictors having a direct effect on Kesi
//   int NV; // Number of Survey Versions
//   int LP; // NUmber of Periods. 
//   array[N_obs] int obs_ind; //  A list of indexes for observed data
//   array[N_est] int est_ind; //  A list of indexes for data where GDP is to be estimated.
//   array[N_eco] int eco_ind; //  A list of indexes for data where GDP is to be estimated.
//   array[M, N_env] int env_ind; // Indexed value for environmetnal predictors of economy
//   array[M] matrix[NP, NP] DmatEnv; // Array of Distance Matrixes for Env Predictors
//   vector[N_env] scale_gdp;
//   array[N_eco] int sver; // GDP per sector
//   // array[N_obs + N_est] int p; // Indicator for particular period
//   array[N_env] int period; // Indicator for particular month
//   matrix[N_eco, K] y; // GDP per sector
//   // For Kesi Predicitions
//   int<lower=0> L; //Number of Years 
//   array[N_env] int ramadan; //  Ramadan
//   array[LP] real P1; // total number of periods
//   matrix[LP, K] period_means;
//   vector[K] k_means;
//   real softplus_alpha;
//   array[DM, N_env] int denv_ind; // Indexed value for environmetnal predictors of Kesi
//   array[DM] matrix[NP, NP] DmatDEnv; // Array of Distance Matrixes for Env Predictors
//   vector[N_env] timber_prices;
//   matrix[L, L] DmatX; // Year distance Matrix
//   array[N_env] int year; //  Time period ID 1:N_est
//   array[N_cases] int kesi; //Number of observed kesi in the time period
// 
// }
// parameters {
//   matrix<lower=0>[N_est, K] gdp_impute;
//   matrix<lower=0>[N_eco, K] gdp_true;
//   array[K] vector[LP] kbp_hat;
// 
//   // Periood Effect
//   array[K] real<lower=0> length_scale_b; // COVAR parameter for periodic GP in Gamma likelihood
//   array[K] real<lower=0> scalesq_b; // Scale parameter for periodic GP in Gamma likelihood
//   // PERIOD SPECIFIC ENVIRONMENTAL EFFECTS ON KESI 
//   array[M, K] real<lower=0> etasq_EnvP;
//   array[M, K] real<lower=0> rhosq_EnvP;
//   array[M, K] real<lower=0> sigma_EnvP;
//   array[LP,M, K] vector[NP] zEnvP;
//   array[M, K] real<lower=0>  cov_EnvP;
//   array[M, K] vector[NP] mu_EnvP;
//   array[K] real<lower=0> sigma_r;
//   array[K] vector[2] br;
//   array[K] real<lower=0> phi;      // Par2 for Gamma
//   matrix[NV, K]   z_ver;     // raw random effect for household survey Version 
//   cholesky_factor_corr[NV] L_ver; // Cholesky factor for Survey Version
//   vector<lower=0>[NV] Sigma_ver;
//   vector[K] gamma_mu;
//   // Kesi Prediction
//   real a;
//   real btimber;
//   vector[K] bgdp;
//   real<lower=0> sigma_r_k;
//   vector[2] br_k;
//   real<lower=0> sigmaX;   // var parameter for Periodic exp_sq in Poisson Likelihood
//   real<lower=0> etasq;   // Scale parameter for Periodic exp_sq in Poisson Likelihood
//   vector[NV] b_sver_kesi;
//   real<lower = 0> sigma_sver_kesi;
//   real<lower=0> omega;
//   real<lower=0> sigma_gdp;
//   real mu_gdp;
//   array[DM] real<lower=0> etasq_DEnv;
//   array[DM] real<lower=0> rhosq_DEnv;
//   array[DM] real<lower=0> sigma_DEnv;
//   array[DM] vector[NP] zDEnv;
//   matrix[K, N_eco] zP;
//   real<lower=0> rhosq;   // COVAR parameter for Periodic exp_sq in Poisson Likelihood
//   real<lower=0> length_scale; // COVAR parameter for Periodic GP in Poisson Likelihood
//   real<lower=0> scalesq;  // Scale parameter for Periodic GP in Poisson Likelihood
//   vector[L] zX;
//   vector[LP] kp_hat;
//   real target_var;
// 
// 
// }
// transformed parameters {
//   matrix[N_est, K] gdp_impute_mu;
//   matrix[N_eco, K] gdp_error;
//   matrix[K, NV]   b_ver;  
//   // Periodic Effect of Sectors
//   array[K] vector[LP] kbp;
//   array[K] cov_matrix[LP] SIGMABP;
//   // Period Specific Effect of ENV on Sectors
//   array[LP, M, K] vector[NP] kEnvP;
//   array[M, K] matrix[NP, NP] LEnv_sigmaP;
//   array[M, K] matrix[NP, NP] SIGMAEnvP;
//    // Direct Effect on Kesi 
//   array[DM] vector[NP] kDEnv;          
//   array[DM] matrix[NP, NP] LDEnv_sigma;
//   array[DM] matrix[NP, NP] SIGMADEnv;
//   cov_matrix[LP] SIGMAP;
//   vector[L] kX;
//   matrix[L, L] LX_sigma;
//   matrix[L, L] SIGMAX;
//   vector[LP] kp;
//   // Environment on Sector Income by period
//   for (m in 1:M){
//     for (k in 1 : K){
//       for(p in 1:LP) {
//         SIGMAEnvP[m,k] = cov_GPL2(DmatEnv[m], etasq_EnvP[m, k], rhosq_EnvP[m, k],
//                                   sigma_EnvP[m, k]);
//         LEnv_sigmaP[m, k] = cholesky_decompose(SIGMAEnvP[m, k]);
//         kEnvP[p,m,k] = LEnv_sigmaP[m,k] * (zEnvP[p,m,k]*cov_EnvP[m, k] + mu_EnvP[m, k]); // This is where the magic happens
//       }
//     }
//   }
//   b_ver = (diag_pre_multiply(Sigma_ver, L_ver) * z_ver)';
//   // Periodic effects for Sectors
//   for (i in 1 : K) {
//     SIGMABP[i] = gp_periodic_cov(P1, scalesq_b[i], length_scale_b[i], 26.0);
//     kbp[i] = cholesky_decompose(SIGMABP[i]) * kbp_hat[i];
//   }
//   // Imputate gdp values based on environmental predictors and training data parameters for all previous periods from 2011 to 2022
//   for (i in 1 : N_est) {
//     for (k in 1 : K) {
//       vector[M] environment_effectP;
//       for(m in 1:M){
//         environment_effectP[m] = kEnvP[period[est_ind[i]], m, k][env_ind[m,est_ind[i]]];
//       }
//         gdp_impute_mu[i, k] = softplus(gamma_mu[k]
//                                   + (kbp[k][period[est_ind[i]]]
//                                   + sum(environment_effectP)
//                                   + br[k][ramadan[est_ind[i]]] * sigma_r[k])
//                                   *scale_gdp[est_ind[i]], softplus_alpha);
//     }
//   }
//   // For the dealing with erros in survey values  
//   for (i in 1 : N_eco) {
//     for (k in 1 : K) {
//       gdp_error[i, k] = softplus(y[i,k] - b_ver[k,sver[i]], softplus_alpha);
//     }
//   }
//   // Periodic effects for Kesi
//   SIGMAP = gp_periodic_cov(P1, scalesq, length_scale, 26.0);
//   kp = cholesky_decompose(SIGMAP) * kp_hat;
//   // Compose GPs
//   SIGMAX = cov_GPL2(DmatX, etasq, rhosq, sigmaX);
//   LX_sigma = cholesky_decompose(SIGMAX);
//   kX = LX_sigma * zX;
//   // Direct effect of Environment on Cases
//   for(m in 1:DM){
//     SIGMADEnv[m] = cov_GPL2(DmatDEnv[m], etasq_DEnv[m],
//                               rhosq_DEnv[m], sigma_DEnv[m]);
//     LDEnv_sigma[m] = cholesky_decompose(SIGMADEnv[m]);
//     kDEnv[m] = LDEnv_sigma[m] * zDEnv[m];
//   }
//   
// }
// model {
//   matrix[N_eco, K] mu;
//   matrix[N_env, K] gdp_merged;
//   vector[N_cases] lambda;
//   to_vector(length_scale_b) ~ normal(0, .2);
//   to_vector(scalesq_b) ~ normal(1, 1);
//   for(k in 1:K)   gamma_mu[k] ~ normal(k_means[k], 1);
//   for(i in 1:LP){
//      for(k in 1:K){
//        //kbp_hat[k][i] ~ normal(period_means[i,k], 1);  
//        kbp_hat[k][i] ~ normal(0, 2);  
//      }
//    }
//   // Ramadan
//   for (i in 1 : K) {
//     br[i] ~ normal(0, 3);
//     sigma_r[i] ~ exponential(1);
//   }
//   // Gamma Scale
//   for(i in 1:K) phi[i] ~ normal(mean(y[,i])/variance(y[,i])+1.5, .3);
//   // Version
//   L_ver ~  lkj_corr_cholesky(2);
//   Sigma_ver ~ exponential(1);
//   to_vector(z_ver) ~ normal(0, 3);
//   // Period Specific Environmental Effects of Kesi
//   for (m in 1:M){
//     for (k in 1 : K){
//       rhosq_EnvP[m,k] ~ exponential(1);
//       etasq_EnvP[m,k] ~ exponential(1);
//       sigma_EnvP[m,k] ~ exponential(1);
//       cov_EnvP[m,k] ~ exponential(1);
//       mu_EnvP[m,k] ~ normal(0, 1);
//       for(p in 1:LP){
//         to_vector(zEnvP[p,m,k]) ~ normal(0, 1);
//       }
//     }
//   }
//   
//   // KESI
//   to_vector(kp_hat) ~ normal(0, 1);
//   to_vector(zP) ~ normal(0, 1);
//   rhosq ~ exponential(1);
//   etasq ~ exponential(1);
//   sigmaX  ~ exponential(1);
//   // Periodic Gausssian Process 
//   length_scale ~ normal(0, .1);
//   scalesq ~ normal(1, 1);
//   to_vector(zX) ~ normal(0, 1);
//   // Environmental GP
//   for(m in 1:DM){
//     rhosq_DEnv[m] ~ exponential(1);
//     etasq_DEnv[m] ~ exponential(1);
//     sigma_DEnv[m] ~ exponential(1);
//     to_vector(zDEnv[m]) ~ normal(0, 1);
//   }
//   
//   br_k ~ normal(0, 1);
//   sigma_r_k ~ exponential(1);
//   btimber ~ normal(0, 1);
//   L_ver ~  lkj_corr_cholesky(2);
//   Sigma_ver ~ exponential(1);
//   to_vector(z_ver) ~ normal(0, 3);
//   b_sver_kesi ~ normal(0, 1);
//   sigma_sver_kesi ~ exponential(1);
//   a ~ normal(12, 3);
//   omega ~ exponential(1);
//   sigma_gdp ~ exponential(1);
//   mu_gdp ~ normal(0, 10);
//   for(i in 1:K) bgdp[i] ~ normal(0, 10);
// 
// 
//   
//   
//   // Estimate parameters for the impact of environmental variabales on each economic sector
//   for (i in 1 : N_eco) {
//     for (k in 1 : K) {
//       vector[M] environment_effectP;
//       for(m in 1:M){
//         environment_effectP[m] = kEnvP[period[eco_ind[i]], m, k][env_ind[m,eco_ind[i]]];
//       }
//       mu[i, k] = softplus(gamma_mu[k] 
//                      + kbp[k][period[eco_ind[i]]] 
//                      + sum(environment_effectP)
//                      + br[k][ramadan[eco_ind[i]]] * sigma_r[k] + b_ver[k, sver[i]],
//                      softplus_alpha);
//     }
//   }
//   for (i in 1 : K) {
//     vector[N_eco] alpha;
//     alpha = mu[ : , i] * phi[i];
//     y[ : , i] ~ gamma(alpha, phi[i]);
//   }
//   for (i in 1 : K) {
//     gdp_impute[ : , i] ~ gamma(gdp_impute_mu[ : , i] * phi[i], phi[i]);
//   }
//   
//   // Deal with measurement error caused by suvey verisons
//     for(t in 1:N_eco){
//     for (k in 1 : K) {
//       real wanted_var =sqrt((y[t,k]*target_var)/y[t,k]); 
// 
//       gdp_true[ t , k] ~ gamma(gdp_error[t ,k] * wanted_var, wanted_var);
//     }
//   }
//   
//   // Assign Imputed Values to Merged
//   for(i in 1:N_est){
//     for(k in 1:K){
//           gdp_merged[i,k] = gdp_impute[i,k];
//      }
//   }
//   // Assign Real Values to Merged
//   for(i in 1:N_eco){
//     for(k in 1:K){
//       gdp_merged[eco_ind[i],k]=gdp_true[ i , k];
//     }
//   }
// 
//  for (i in 1 : N_cases) {
//     vector[DM] direct_environmental_effect;
//     vector[K]  economic_impact;
//     for(m in 1:DM){
//               direct_environmental_effect[m] = kDEnv[m][denv_ind[m, i]];
//     }
//      for(k in 1:K){
//       economic_impact[k] = (bgdp[k]*sigma_gdp + mu_gdp)*log1p(gdp_merged[i,k]/10);
//      }
//     lambda[i] =     softplus(a +
//                     kp[period[i]] +
//                     sum(economic_impact) +
//                     kX[year[i]] +
//                     br_k[ramadan[i]] * sigma_r_k +
//                     sum(direct_environmental_effect) +
//                     btimber*timber_prices[i], softplus_alpha);
//   }
//   // Estimate the Number of Kesi!  
//   for (i in 1 : N_cases) {
//     kesi[i] ~ neg_binomial_2(lambda[i], omega);
//   }
// }
// 
// 
