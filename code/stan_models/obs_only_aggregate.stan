// file: kesi_direct_from_y_with_ME.stan
functions {
  // Smoothed Exponential gaussian Process (unchanged)
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
  // Period Gaussian Process (unchanged)
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
  // Stable softplus (unchanged)
  real softplus(real x, real alpha) {
    return log1p(exp(x * alpha)) / alpha;
  }
}

data {
  // Observations aligned with y
  int<lower=1> N_cases;             // rows of y
  int<lower=1> K;                   // sectors (cols of y)
  matrix[N_cases, K] y;             // observed GDP per sector (as provided)
  array[N_cases] int<lower=0> kesi; // case counts
  vector[N_cases] timber_prices;

  // Time indexing
  int<lower=1> L;                   // years
  array[N_cases] int<lower=1,upper=L>  year;
  int<lower=1> LP;                  // periods
  array[LP] real P1;
  real target_var;
  array[N_cases] int<lower=1,upper=LP> period;
  array[N_cases] int<lower=1,upper=2>  ramadan;

  // Year GP distance
  matrix[L, L] DmatX;

  // Direct environmental effects on kesi
  int<lower=0> DM;
  int<lower=0> NP;
  array[DM] matrix[NP, NP] DmatDEnv;
  array[DM, N_cases] int<lower=1,upper=NP> denv_ind;

  // Survey version info for y
  int<lower=1> NV;
  array[N_cases] int<lower=1,upper=NV> sver;

  // Link stability
  real<lower=0> softplus_alpha;
}

parameters {
  // Periodic GP for kesi
  vector[LP] kp_hat;
  real<lower=0> rhosq;
  real<lower=0> length_scale;
  real<lower=0> scalesq;

  // Year GP
  vector[L] zX;
  real<lower=0> sigmaX;
  real<lower=0> etasq;

  // Direct env GPs
  array[DM] real<lower=0> etasq_DEnv;
  array[DM] real<lower=0> rhosq_DEnv;
  array[DM] real<lower=0> sigma_DEnv;
  array[DM] vector[NP] zDEnv;

  // Survey version × sector adjustment on y
  matrix[NV, K]   z_ver;            // raw effects
  cholesky_factor_corr[NV] L_ver;   // version correlation
  vector<lower=0>[NV] Sigma_ver;    // version scales

  // --- Measurement error latent true GDP (positive) ---
  matrix<lower=0>[N_cases, K] y_true;   // resampled latent GDP after version adjustment
  //array[K] real<lower=0> tau_y;                  // precision for ME layer (higher -> tighter around adj y)

  // Fixed effects for kesi
  real a;
  real btimber;
  real bgdp;
  vector[2] br_k;

  // NB2 dispersion
  real<lower=0> omega;
}

transformed parameters {
  // Year GP
  matrix[L, L] SIGMAX = cov_GPL2(DmatX, etasq, rhosq, sigmaX);
  matrix[L, L] LX_sigma = cholesky_decompose(SIGMAX);
  vector[L] kX = LX_sigma * zX;

  // Direct env GPs
  array[DM] vector[NP] kDEnv;
  for (m in 1:DM) {
    matrix[NP, NP] S = cov_GPL2(DmatDEnv[m], etasq_DEnv[m], rhosq_DEnv[m], sigma_DEnv[m]);
    kDEnv[m] = cholesky_decompose(S) * zDEnv[m];
  }

  // Survey version adjustment on y (K × NV), same as your working code
  matrix[K, NV] b_ver = (diag_pre_multiply(Sigma_ver, L_ver) * z_ver)';

  // Periodic GP for kesi
  cov_matrix[LP] SIGMAP = gp_periodic_cov(P1, scalesq, length_scale, 26.0);
  vector[LP] kp = cholesky_decompose(SIGMAP) * kp_hat;
}

model {
  // Priors (unchanged scales)
  to_vector(kp_hat) ~ normal(0, 1);
  rhosq ~ exponential(1);
  etasq ~ exponential(1);
  sigmaX ~ exponential(1);
  length_scale ~ normal(0, .1);
  scalesq ~ normal(1, .5);
  to_vector(zX) ~ normal(0, 1);

  for (m in 1:DM) {
    rhosq_DEnv[m] ~ exponential(1);
    etasq_DEnv[m] ~ exponential(1);
    sigma_DEnv[m] ~ exponential(1);
    to_vector(zDEnv[m]) ~ normal(0, 1);
  }

  // Survey version hierarchy
  L_ver ~ lkj_corr_cholesky(2);
  Sigma_ver ~ exponential(1);
  to_vector(z_ver) ~ normal(0, 3);

  // Kesi fixed effects
  bgdp ~ normal(0, 1);
  btimber ~ normal(0, 1);
  br_k ~ normal(0, 1);
  a ~ normal(28, 3);
  omega ~ exponential(1);

  // --- Measurement error layer on version-adjusted y ---
  // Define version-adjusted mean (positive via softplus), then draw y_true ~ Gamma(mean=adj_y, var=adj_y/tau_y)
  //tau_y ~ normal(10, 1);  // you can swap to target_var-driven if you prefer
  for (i in 1:N_cases) {
    for (k in 1:K) {
      real y_adj = softplus( y[i, k] - b_ver[k, sver[i]], softplus_alpha );  // version-corrected, positive
      // Gamma with mean = y_adj, variance = y_adj / tau_y -> shape = y_adj * tau_y, rate = tau_y
      real wanted_var =sqrt((y[i,k]*target_var)/y[i,k]); 
      y_true[i, k] ~ gamma( y_adj * wanted_var, wanted_var );
    }
  }
  
  
  // --- Kesi likelihood using log(y_true) ---
  for (i in 1:N_cases) {
    vector[DM] dE;
    for (m in 1:DM) dE[m] = kDEnv[m][ denv_ind[m, i] ];
    // log of positive latent y_true
    real log_y_i = ( log(sum( y[i,] ) ));
    real f = a
             + kp[ period[i] ]
             + (log_y_i * bgdp)
             + kX[ year[i] ]
             + br_k[ ramadan[i] ]
             + sum(dE)
             + btimber * timber_prices[i];

    kesi[i] ~ neg_binomial_2( softplus(f, softplus_alpha), omega );
  }
}
// 
// generated quantities {
//   vector[N_cases] log_lik;
//   vector[N_cases] lambda;
// 
//   for (i in 1:N_cases) {
//     // rebuild linear predictor using already-sampled y_true
//     row_vector[K] log_y_i = log( to_row_vector( y_true[i] ) );
// 
//     // periodic piece rederived (cheap & deterministic given params)
//     cov_matrix[LP] SIGMAP = gp_periodic_cov(P1, scalesq, length_scale, 26.0);
//     vector[LP] kp = cholesky_decompose(SIGMAP) * kp_hat;
// 
//     real f = a
//              + kp[ period[i] ]
//              + (log_y_i * bgdp)
//              + kX[ year[i] ]
//              + br_k[ ramadan[i] ]
//              + btimber * timber_prices[i];
// 
//     // add direct env effects
//     {
//       real add_env = 0;
//       for (m in 1:DM) add_env += (cholesky_decompose(cov_GPL2(DmatDEnv[m], etasq_DEnv[m], rhosq_DEnv[m], sigma_DEnv[m])) * zDEnv[m])[ denv_ind[m, i] ];
//       f += add_env;
//     }
// 
//     lambda[i] = softplus(f, softplus_alpha);
//     log_lik[i] = neg_binomial_2_lpmf(kesi[i] | lambda[i], omega);
//   }
// }
