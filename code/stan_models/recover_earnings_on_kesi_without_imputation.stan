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

data{

  // For GDP Estimation
  int<lower=0> N_cases;  // Number of Observed Cases for Prediction
  int<lower=0> N_obs;  // Number of time points where cases are observed and we have comparable income data
  int<lower=0> N_est;  // Number of time points GDP data is unavailable but cases are
  int<lower=0> N_eco;  // Number of economic observations
  int<lower=0> N_env;  // Number of environmental observations
  int<lower=0> K;  // Number of sectors
  int<lower=0> NP;  // Number of Percentiles for GPs
  matrix[N_cases,K] y2;
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
     
parameters{
     real a;
     vector[K] bgdp;
     real mu_gdp;
     real<lower = 0> sigma_gdp;
     real<lower = 0> omega;
     vector[L] zX;  // Latent variables for year effect
     vector[LP] kp_hat;  // Raw periodic effects for Kesi
     matrix[K, N_eco] zP;  // Latent variables for periodic GP
     real<lower=0> rhosq;  // Covariance parameter for periodic GP in Poisson likelihood
     real<lower=0> length_scale;  // Length scale for periodic GP in Poisson likelihood
     real<lower=0> scalesq;  // Scale parameter for periodic GP in Poisson likelihood
     array[DM] real<lower=0> etasq_DEnv;  // Scale parameters for direct environmental effects
     array[DM] real<lower=0> rhosq_DEnv;  // Length scales for direct environmental effects
     array[DM] real<lower=0> sigma_DEnv;  // Noise terms for direct environmental effects
     array[DM] vector[NP] zDEnv;  // Latent variables for direct environmental effects
     real btimber;  // Timber price effect
     real<lower=0> sigma_r_k;  // Standard deviation for Kesi Ramadan effect
     vector[2] br_k;  // Ramadan effect parameters for Kesi
     real<lower=0> sigmaX;  // Variance parameter for periodic GP in Poisson likelihood
     real<lower=0> etasq;  // Scale parameter for periodic GP in Poisson likelihood
}

transformed parameters {
  // Direct Effect on Kesi 
  array[DM] vector[NP] kDEnv;  // Direct environmental effects
  array[DM] matrix[NP, NP] LDEnv_sigma;  // Cholesky factors for direct environmental effects
  array[DM] matrix[NP, NP] SIGMADEnv;  // Covariance matrices for direct environmental effects
  
  cov_matrix[LP] SIGMAP;  // Covariance matrix for periodic GP
  vector[L] kX;  // Year effect
  matrix[L, L] LX_sigma;  // Cholesky factor for year effect
  matrix[L, L] SIGMAX;  // Covariance matrix for year effect
  vector[LP] kp;  // Periodic effects for Kesi
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

model{
    vector[N_cases] lambda;
    bgdp ~ normal( 0 , 5 );
    mu_gdp ~ normal( 0 , 5 );
    sigma_gdp ~ exponential(1);
    //a ~ normal(mean(kesi)+max(y2), 1);
    a ~ normal(mean(kesi), 1);
    omega ~ exponential(1);
    // Priors for periodic GP parameters
    length_scale ~ normal(0, .1);
    scalesq ~ normal(1, 1);
    to_vector(zX) ~ normal(0, 1);

    to_vector(kp_hat) ~ normal(0, 1);
    to_vector(zP) ~ normal(0, 1);
    rhosq ~ exponential(1);
    etasq ~ exponential(1);
    sigmaX ~ exponential(1);
    
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

    for ( i in 1:N_cases ) {
      vector[DM] direct_environmental_effect;
      for (m in 1:DM) {
        direct_environmental_effect[m] = kDEnv[m][denv_ind[m, i]];
      }
        lambda[i] = softplus(a +  kp[period[i]] + kX[year[i]] +
                    br_k[ramadan[i]] * sigma_r_k +
                    sum(direct_environmental_effect) +
                    btimber * timber_prices[i] +
                    dot_product((bgdp * sigma_gdp + mu_gdp), log1p(to_vector(y2[i,]))), softplus_alpha);
        
    }
      // Likelihood for observed Kesi
    for (i in 1:N_cases) {
      kesi[i] ~ neg_binomial_2(lambda[i], omega);
    }
}

