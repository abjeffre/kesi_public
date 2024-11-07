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
  int DM; // Number of Environmental Predictors having a direct effect on Kesi
  int NM; // NUmber of Months
  int NV; // Number of Survey Versions
  int LP; // NUmber of Periods. 
  array[N_obs] int obs_ind; //  A list of indexes for observed data
  array[N_est] int est_ind; //  A list of indexes for data where GDP is to be estimated.
  array[N_eco] int eco_ind; //  A list of indexes for data where GDP is to be estimated.
  array[M, N_env] int env_ind; // Indexed value for environmetnal predictors of economy
  array[DM, N_env] int denv_ind; // Indexed value for environmetnal predictors of Kesi
  array[M] matrix[NP, NP] DmatEnv; // Array of Distance Matrixes for Env Predictors
  array[DM] matrix[NP, NP] DmatDEnv; // Array of Distance Matrixes for Env Predictors
  vector[N_env] scale_gdp;
  vector[N_env] timber_prices;
  array[N_eco] int sver; // GDP per sector

  // array[N_obs + N_est] int p; // Indicator for particular period
  array[N_env] int period; // Indicator for particular month
  matrix[N_eco, K] y; // GDP per sector
  // For Kesi Predicitions
  int<lower=0> L; //Number of Years 
  array[N_cases] int kesi; //Number of observed kesi in the time period
  array[N_env] int year; //  Time period ID 1:N_est
  array[N_env] int ramadan; //  Ramadan
  matrix[L, L] DmatX; // Year distance Matrix
  array[LP] real P1; // total number of periods
  vector[N_eco] gdp; // observed GDP 
  matrix[LP, K] period_means;
  real softplus_alpha;
  real target_var;
}
parameters {
  matrix<lower=0>[N_est, K] gdp_impute;
  matrix<lower=0>[N_eco, K] gdp_true;
  vector[LP] kp_hat;
  array[K] vector[LP] kbp_hat;
  vector[L] zX;

  // Periood Effect
  matrix[K, N_eco] zP;
  real<lower=0> rhosq;   // COVAR parameter for Periodic exp_sq in Poisson Likelihood
  real<lower=0> length_scale; // COVAR parameter for Periodic GP in Poisson Likelihood
  real<lower=0> scalesq;  // Scale parameter for Periodic GP in Poisson Likelihood
  array[K] real<lower=0> length_scale_b; // COVAR parameter for periodic GP in Gamma likelihood
  array[K] real<lower=0> scalesq_b; // Scale parameter for periodic GP in Gamma likelihood

  // PERIOD SPECIFIC ENVIRONMENTAL EFFECTS ON KESI 
  array[M, K] real<lower=0> etasq_EnvP;
  array[M, K] real<lower=0> rhosq_EnvP;
  array[M, K] real<lower=0> sigma_EnvP;
  array[LP,M, K] vector[NP] zEnvP;
  array[M, K] real<lower=0>  cov_EnvP;
  array[M, K] vector[NP] mu_EnvP;

  // DIRECT EFFECTS ON KASES
  array[DM] real<lower=0> etasq_DEnv;
  array[DM] real<lower=0> rhosq_DEnv;
  array[DM] real<lower=0> sigma_DEnv;
  array[DM] vector[NP] zDEnv;
  // Kesi Prediction
  real a;
  real btimber;
  vector[K] bgdp;
  array[K] real<lower=0> sigma_r;
  array[K] vector[2] br;
  real<lower=0> sigma_r_k;
  vector[2] br_k;
  real<lower=0> sigmaX;   // var parameter for Periodic exp_sq in Poisson Likelihood
  real<lower=0> etasq;   // Scale parameter for Periodic exp_sq in Poisson Likelihood
  array[K] real<lower=0> phi;      // Par2 for Gamma
  matrix[NV, K]   z_ver;     // raw random effect for household survey Version 
  cholesky_factor_corr[NV] L_ver; // Cholesky factor for Survey Version
  vector<lower=0>[NV] Sigma_ver;
  vector[NV] b_sver_kesi;
  real<lower = 0> sigma_sver_kesi;
  real<lower=0> omega;
  real<lower=0> sigma_gdp;
  real mu_gdp;

}
transformed parameters {
  matrix[N_est, K] gdp_impute_mu;
  matrix[N_eco, K] gdp_error;
  cov_matrix[LP] SIGMAP;
  vector[L] kX;
  matrix[L, L] LX_sigma;
  matrix[L, L] SIGMAX;
  vector[LP] kp;
  matrix[K, NV]   b_ver;  

  // Periodic Effect of Sectors
  array[K] vector[LP] kbp;
  array[K] cov_matrix[LP] SIGMABP;
  // Period Specific Effect of ENV on Sectors
  array[LP, M, K] vector[NP] kEnvP;
  array[M, K] matrix[NP, NP] LEnv_sigmaP;
  array[M, K] matrix[NP, NP] SIGMAEnvP;
 // Direct Effect on Kesi 
  array[DM] vector[NP] kDEnv;          
  array[DM] matrix[NP, NP] LDEnv_sigma;
  array[DM] matrix[NP, NP] SIGMADEnv;
  // Compose GPs
  SIGMAX = cov_GPL2(DmatX, etasq, rhosq, sigmaX);
  LX_sigma = cholesky_decompose(SIGMAX);
  kX = LX_sigma * zX;
  // Direct effect of Environment on Cases
  for(m in 1:DM){
    SIGMADEnv[m] = cov_GPL2(DmatDEnv[m], etasq_DEnv[m],
                              rhosq_DEnv[m], sigma_DEnv[m]);
    LDEnv_sigma[m] = cholesky_decompose(SIGMADEnv[m]);
    kDEnv[m] = LDEnv_sigma[m] * zDEnv[m];
  }
  // Environment on Sector Income by period
  for (m in 1:M){
    for (k in 1 : K){
      for(p in 1:LP) {
        SIGMAEnvP[m,k] = cov_GPL2(DmatEnv[m], etasq_EnvP[m, k], rhosq_EnvP[m, k],
                                  sigma_EnvP[m, k]);
        LEnv_sigmaP[m, k] = cholesky_decompose(SIGMAEnvP[m, k]);
        kEnvP[p,m,k] = LEnv_sigmaP[m,k] * (zEnvP[p,m,k,]*cov_EnvP[m, k] + mu_EnvP[m, k,]); // This is where the magic happens
      }
    }
  }

  b_ver = (diag_pre_multiply(Sigma_ver, L_ver) * z_ver)';
  // Periodic effects for Kesi
  SIGMAP = gp_periodic_cov(P1, scalesq, length_scale, 26.0);
  kp = cholesky_decompose(SIGMAP) * kp_hat;
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
        gdp_impute_mu[i, k] = softplus( 
                                  (kbp[k][period[est_ind[i]]]
                                  + sum(environment_effectP)
                                  + br[k][ramadan[est_ind[i]]])
                                  , softplus_alpha);
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
  vector[N_cases] lambda;
  to_vector(kp_hat) ~ normal(0, 1);
  rhosq ~ exponential(1);
  etasq ~ exponential(1);
  sigmaX  ~ exponential(1);
  // Periodic Gausssian Process 
  length_scale ~ normal(0, .1);
  scalesq ~ normal(1, .5);
  for(i in 1:LP){
     for(k in 1:K){
       kbp_hat[k][i] ~ normal(period_means[i,k], 1);    
     }
   }
  to_vector(length_scale_b) ~ normal(0, .1);
  to_vector(scalesq_b) ~ normal(1, .5);
  to_vector(zX) ~ normal(0, 1);
  // Environmental GP
  for(m in 1:DM){
    rhosq_DEnv[m] ~ exponential(1);
    etasq_DEnv[m] ~ exponential(1);
    sigma_DEnv[m] ~ exponential(1);
    to_vector(zDEnv[m]) ~ normal(0, 1);
  }
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
  // Ramadan
  for (i in 1 : K) {
    br[i] ~ normal(0, 1);
    sigma_r[i] ~ exponential(1);
  }
  br_k ~ normal(0, 1);
  sigma_r_k ~ exponential(1);
  // Income effect
  for(i in 1:K) bgdp[i] ~ normal(0, 2);
  // Gamma Scale
  for(i in 1:K) phi[i] ~ normal(mean(y[,i])/variance(y[,i])+1.5, .1);
  // Timber
  btimber ~ normal(0, 1);
  L_ver ~  lkj_corr_cholesky(2);
  Sigma_ver ~ exponential(1);
  to_vector(z_ver) ~ normal(0, 3);
  b_sver_kesi ~ normal(0, 1);
  sigma_sver_kesi ~ exponential(1);
  a ~ normal(16, .5); // Note this prior is chosen given testing using frequentist models to determine likely range. 
  omega ~ exponential(1);
  sigma_gdp ~ exponential(1);
  mu_gdp ~ normal(0, 1);

  // Estimate parameters for the impact of environmental variabales on each economic sector
  for (i in 1 : N_eco) {
    for (k in 1 : K) {
      vector[M] environment_effectP;
      for(m in 1:M){
        environment_effectP[m] = kEnvP[period[eco_ind[i]], m, k][env_ind[m,eco_ind[i]]];
      }
      mu[i, k] = softplus( 
                       kbp[k][period[eco_ind[i]]] 
                     + sum(environment_effectP)
                     + br[k][ramadan[eco_ind[i]]] 
                     + b_ver[k, sver[i]],
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
      real wanted_var =sqrt((y[t,k]*target_var)/y[t,k]); 

      gdp_true[ t , k] ~ gamma(gdp_error[t ,k] * wanted_var, wanted_var);
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
  for (i in 1 : N_cases) {
    vector[DM] direct_environmental_effect;
    vector[K]  economic_impact;
    for(m in 1:DM){
              direct_environmental_effect[m] = kDEnv[m][denv_ind[m, i]];
    }
    lambda[i] =     softplus(a +
                    kp[period[i]] +
                    dot_product((bgdp*sigma_gdp + mu_gdp), to_vector(log1p(gdp_merged[i,]))) +
                    kX[year[i]] +
                    br_k[ramadan[i]] +
                    sum(direct_environmental_effect) +
                    btimber*timber_prices[i], softplus_alpha);
  }
  // Estimate the Number of Kesi!  
  for (i in 1 : N_cases) {
    kesi[i] ~ neg_binomial_2(lambda[i], omega);
  }
}

generated quantities{
  vector[N_cases] log_lik;
  vector[N_cases] lambda;
  matrix[N_env, K] gdp_merged;

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

  for (i in 1 : N_cases) {
    vector[DM] direct_environmental_effect;
    vector[K]  economic_impact;
    for(m in 1:DM){
              direct_environmental_effect[m] = kDEnv[m][denv_ind[m, i]];
    }
    lambda[i] =     softplus(a +
                    kp[period[i]] +
                    dot_product((bgdp*sigma_gdp + mu_gdp), to_vector(log1p(gdp_merged[i,]))) +
                    kX[year[i]] +
                    br_k[ramadan[i]] +
                    sum(direct_environmental_effect) +
                    btimber*timber_prices[i], softplus_alpha);
  }

 for ( i in 1:N_cases ) {
  log_lik[i] = neg_binomial_2_lpmf( kesi[i] | lambda[i], omega );
  }
}
