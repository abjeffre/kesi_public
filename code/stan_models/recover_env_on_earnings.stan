functions {
  // Smoothed Exponential gaussian Process
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1 : (N - 1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1) : N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i, j]));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
  // Period Gaussian Process
  matrix cov_periodic(matrix x, real sq_alpha, real sq_rho, real delta,
                      real period) {
    int NOBS = dims(x)[1];
    matrix[NOBS, NOBS] K;
    for (i in 1 : (NOBS - 1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1) : NOBS) {
        K[i, j] = (sq_alpha ^ 2)
                  * exp(-2 * (sin(pi() * x[i, j] / period) ^ 2)
                        / (sq_rho ^ 2));
        K[j, i] = K[i, j];
      }
    }
    K[NOBS, NOBS] = sq_alpha + delta;
    return K;
  }
  // Merge missing for imputation 
  real softplus(real x, real alpha) {
    real output;
    output =log1p(exp(x*alpha))/alpha; 
    return output;
  }
}
data {
  // For GDP Estimation
  int<lower=0> N_cases; // Number of Observed Cases for Prediction
  int<lower=0> N_obs; // Number of time where cases are obsereved and we have comparable income data
  int<lower=0> N_est; // Number of time points gdp data is unavailible but cases are.
  int<lower=0> N_eco; // Number of economic obervations
  int<lower=0> N_env; // Number of economic obervations
  int<lower=0> K; // Number of sectors
  int<lower=0> NP; // Number of Percentiles for GPs. 
  int M; // Number of Environmental Predictors
  int NV; // Number of Survey Versions
  int LP; // NUmber of Periods. 
  array[N_obs] int obs_ind; //  A list of indexes for observed data
  array[N_est] int est_ind; //  A list of indexes for data where GDP is to be estimated.
  array[N_eco] int eco_ind; //  A list of indexes for data where GDP is to be estimated.
  array[M, N_env] int env_ind; // Indexed value for environmetnal predictors of economy
  array[M] matrix[NP, NP] DmatEnv; // Array of Distance Matrixes for Env Predictors
  vector[N_env] scale_gdp;
  array[N_eco] int sver; // GDP per sector
  // array[N_obs + N_est] int p; // Indicator for particular period
  array[N_env] int period; // Indicator for particular month
  matrix[N_eco, K] y; // GDP per sector
  // For Kesi Predicitions
  int<lower=0> L; //Number of Years 
  array[N_env] int ramadan; //  Ramadan
  array[LP] real P1; // total number of periods
  matrix[LP, K] period_means;
  vector[K] k_means;
  real softplus_alpha;
  real target_var;  // Target variance
}
parameters {
  matrix<lower=0>[N_est, K] gdp_impute;
  matrix<lower=0>[N_eco, K] gdp_true;
  array[K] vector[LP] kbp_hat;

  // Periood Effect
  array[K] real<lower=0> length_scale_b; // COVAR parameter for periodic GP in Gamma likelihood
  array[K] real<lower=0> scalesq_b; // Scale parameter for periodic GP in Gamma likelihood
  // PERIOD SPECIFIC ENVIRONMENTAL EFFECTS ON KESI 
  array[M, K] real<lower=0> etasq_EnvP;
  array[M, K] real<lower=0> rhosq_EnvP;
  array[M, K] real<lower=0> sigma_EnvP;
  array[LP,M, K] vector[NP] zEnvP;
  array[M, K] real<lower=0>  cov_EnvP;
  array[M, K] vector[NP] mu_EnvP;
  // Kesi Prediction
  array[K] real<lower=0> sigma_r;
  array[K] vector[2] br;
  array[K] real<lower=0> phi;      // Par2 for Gamma
  matrix[NV, K]   z_ver;     // raw random effect for household survey Version 
  cholesky_factor_corr[NV] L_ver; // Cholesky factor for Survey Version
  vector<lower=0>[NV] Sigma_ver;
  vector[K] gamma_mu;

}
transformed parameters {
  matrix[N_est, K] gdp_impute_mu;
  matrix[N_eco, K] gdp_error;
  matrix[K, NV]   b_ver;  
  // Periodic Effect of Sectors
  array[K] vector[LP] kbp;
  array[K] cov_matrix[LP] SIGMABP;
  // Period Specific Effect of ENV on Sectors
  array[LP, M, K] vector[NP] kEnvP;
  array[M, K] matrix[NP, NP] LEnv_sigmaP;
  array[M, K] matrix[NP, NP] SIGMAEnvP;
  // Environment on Sector Income by period
  for (m in 1:M){
    for (k in 1 : K){
      for(p in 1:LP) {
        SIGMAEnvP[m,k] = cov_GPL2(DmatEnv[m], etasq_EnvP[m, k], rhosq_EnvP[m, k],
                                  sigma_EnvP[m, k]);
        LEnv_sigmaP[m, k] = cholesky_decompose(SIGMAEnvP[m, k]);
        kEnvP[p,m,k] = LEnv_sigmaP[m,k] * (zEnvP[p,m,k]*cov_EnvP[m, k] + mu_EnvP[m, k]); // This is where the magic happens
      }
    }
  }
  b_ver = (diag_pre_multiply(Sigma_ver, L_ver) * z_ver)';
  // Periodic effects for Sectors
  for (i in 1 : K) {
    SIGMABP[i] = gp_periodic_cov(P1, scalesq_b[i], length_scale_b[i], 26.0);
    kbp[i] = cholesky_decompose(SIGMABP[i]) * kbp_hat[i];
  }
  // Imputate gdp values based on environmental predictors and training data parameters for all previous periods from 2011 to 2022
  for (i in 1 : N_est) {
    for (k in 1 : K) {
      vector[M] environment_effectP;
      for(m in 1:M){
        environment_effectP[m] = kEnvP[period[est_ind[i]], m, k][env_ind[m,est_ind[i]]];
      }
        gdp_impute_mu[i, k] = softplus(gamma_mu[k] +
                                   (kbp[k][period[est_ind[i]]]
                                  + sum(environment_effectP)
                                  + br[k][ramadan[est_ind[i]]] * sigma_r[k])
                                  *scale_gdp[est_ind[i]], softplus_alpha);
    }
  }
  // For the dealing with erros in survey values  
  for (i in 1 : N_eco) {
    for (k in 1 : K) {
      gdp_error[i, k] = softplus(y[i,k] - b_ver[k,sver[i]], softplus_alpha);
    }
  }
}
model {
  matrix[N_eco, K] mu;
  matrix[N_env, K] gdp_merged;
  to_vector(length_scale_b) ~ normal(0, .2);
  to_vector(scalesq_b) ~ normal(1, 1);
  for(k in 1:K)   gamma_mu[k] ~ normal(k_means[k], 1);
  for(i in 1:LP){
     for(k in 1:K){
       kbp_hat[k][i] ~ normal(period_means[i,k], 1);  
       //kbp_hat[k][i] ~ normal(0, 2);  
     }
   }
  // Ramadan
  for (i in 1 : K) {
    br[i] ~ normal(0, 3);
    sigma_r[i] ~ exponential(1);
  }
  // Gamma Scale
  for(i in 1:K) phi[i] ~ normal(mean(y[,i])/variance(y[,i])+1.5, .3);
  // Version
  L_ver ~  lkj_corr_cholesky(2);
  Sigma_ver ~ exponential(1);
  to_vector(z_ver) ~ normal(0, .5);
  // Period Specific Environmental Effects of Kesi
  for (m in 1:M){
    for (k in 1 : K){
      rhosq_EnvP[m,k] ~ exponential(1);
      etasq_EnvP[m,k] ~ exponential(1);
      sigma_EnvP[m,k] ~ exponential(1);
      cov_EnvP[m,k] ~ exponential(1);
      mu_EnvP[m,k] ~ normal(0, 1);
      for(p in 1:LP){
        to_vector(zEnvP[p,m,k]) ~ normal(0, 1);
      }
    }
  }
  // Estimate parameters for the impact of environmental variabales on each economic sector
  for (i in 1 : N_eco) {
    for (k in 1 : K) {
      vector[M] environment_effectP;
      for(m in 1:M){
        environment_effectP[m] = kEnvP[period[eco_ind[i]], m, k][env_ind[m,eco_ind[i]]];
      }
      mu[i, k] = softplus(gamma_mu[k]+ 
                     + kbp[k][period[eco_ind[i]]] 
                     + sum(environment_effectP)
                     + br[k][ramadan[eco_ind[i]]] * sigma_r[k] + b_ver[k, sver[i]],
                     softplus_alpha);
    }
  }
  for (i in 1 : K) {
    vector[N_eco] alpha;
    alpha = mu[ : , i] * phi[i];
    y[ : , i] ~ gamma(alpha, phi[i]);
  }
  for (i in 1 : K) {
    gdp_impute[ : , i] ~ gamma(gdp_impute_mu[ : , i] * phi[i], phi[i]);
  }
  
  // Deal with measurement error caused by suvey verisons
  for(t in 1:N_eco){
    for (k in 1 : K) {
      real wanted_var = sqrt((y[t, k] * target_var) / y[t, k]);
      gdp_true[t, k] ~ gamma(gdp_error[t, k] * wanted_var, wanted_var);

    }
  }
  
  // Assign Imputed Values to Merged
  for(i in 1:N_est){
    for(k in 1:K){
          gdp_merged[i,k] = gdp_impute[i,k];
     }
  }
  // Assign Real Values to Merged
  for(i in 1:N_eco){
    for(k in 1:K){
      gdp_merged[eco_ind[i],k]=gdp_true[ i , k];
    }
  }
}

