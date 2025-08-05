functions {
  /* cumulative-logit log-PDF for a single response
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   sigma: latent sd parameter
   *   thres: ordinal thresholds
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real cumulative_logit_lpmf(int y, real mu, real sigma, vector thres) {
     int nthres = num_elements(thres);
     if (y == 1) {
       return log_inv_logit((thres[1] - mu) / sigma);
     } else if (y == nthres + 1) {
       return log1m_inv_logit((thres[nthres] - mu) / sigma);
     } else {
       return log_diff_exp(
         log_inv_logit((thres[y] - mu) / sigma), 
         log_inv_logit((thres[y - 1] - mu) / sigma)
       );
     }
   }
}
data {
  //basic information and indices
  int<lower=1> I; //number of individuals (#Fband)
  int<lower=1> N; //number of observations
  int<lower=1> nyear; //number of years
  array[N] int<lower=1> id; //index of individual observations
  array[N] int<lower=1> year_id; //index of years
  
  //population-level predictors (non-selection fixed effects)
  int<lower=1> nthres_clutch;
  int<lower=1> nthres_hatch;
  int<lower=1> k_clutch; //number of covariates excluding global intercept
  matrix[N, k_clutch] X_clutch; // population-level design matrix excl. intercept
  int<lower=1> k_hatch;
  matrix[N, k_hatch] X_hatch;
  int<lower=1> k_survive;
  matrix[N, k_survive] X_survive;
  int<lower=1> k_mass;
  matrix[N, k_mass] X_mass;
  
  //environmental measures
  vector[N] dateWSD;
  vector[N] attemptSD;
  
  //response variables
  array[N] int clutch;
  array[N] int hatch;
  array[N] int survive;
  array[N] real mass;
}
transformed data {
  int<lower=1> nresp = 3; //number of response variables (fitness components)
  int<lower=1> ntrait = 2; //number of RN components under selection (intercepts and residuals)
  int ntrait_cor = (ntrait * ntrait - ntrait) / 2; //unique correlational selection effects
  int<lower=1> nselect = ntrait + ntrait + ntrait_cor; //number of selection effects
  //cholesky_factor_corr[nselect] I_mat = cholesky_decompose(diag_matrix(rep_vector(1,nselect)));
  //real delta = -9; //noise added for identification on log scale; exp(-9) = 0.0001
}
parameters {
  //population fixed effects
  ordered[nthres_clutch] mu_c_0; //intercepts
  vector[k_clutch] coef_c; //regression coefficients
  //real v_0;  // disc intercept
  
  ordered[nthres_hatch] mu_h_0; //intercepts 
  vector[k_hatch] coef_h;
  //real<lower=0> shape_h;
  
  real int_c;
  real int_h;
  real mu_s_0;
  vector[k_hatch] coef_s;
  real<lower=0> sd_O_s;
  vector[N] O_s_0;
  
  real mu_m_0;
  vector[k_mass] coef_m;
  real<lower=0> sd_sigma_m;
  
  //selection coefficients
  vector[ntrait] b_h; //linear selection
  vector[ntrait] db_h;
  vector[ntrait + ntrait_cor] q_h; //nonlinear selection
  vector[ntrait + ntrait_cor] dq_h;
  vector[ntrait] b_s;
  vector[ntrait] db_s;
  vector[ntrait + ntrait_cor] q_s;
  vector[ntrait + ntrait_cor] dq_s;
  vector[ntrait] b_m;
  vector[ntrait] db_m;
  vector[ntrait + ntrait_cor] q_m;
  vector[ntrait + ntrait_cor] dq_m;
  
  //random effects for reaction norm
  vector<lower=0>[ntrait] sd_I_c; //RN parameter SDs
  matrix[I, ntrait] std_dev_c; //individual scores (unit scale/z-score deviations)
  cholesky_factor_corr[ntrait] cor_chol; //RN parameter correlations
  
  //random effects for fitness
  //individual variation unexplained by RNs
  real<lower=0> sd_I_h;
  vector[I] std_dev_h;
  real<lower=0> sd_I_s;
  vector[I] std_dev_s;
  real<lower=0> sd_I_m;
  vector[I] std_dev_m;
  
  //random intercepts for year
  real<lower=0> sd_Y_c_0;
  vector[nyear] std_dev_Y_c_0;
  real<lower=0> sd_Y_h_0;
  vector[nyear] std_dev_Y_h_0;
  real<lower=0> sd_Y_s_0;
  vector[nyear] std_dev_Y_s_0;
  real<lower=0> sd_Y_m_0;
  vector[nyear] std_dev_Y_m_0;
}
transformed parameters {
  //random effects are separated into separate z-score, SD, and correlation parameters for 
  //MCMC sampling efficiency. We combine these parameters together to calculate the 
  //VCV matrix from the SDs and corr matrix, and then we use this matrix to scale the
  //individual z-scores, providing random effects on the appropriate scale
  
  matrix[I, ntrait] I_c = std_dev_c * diag_pre_multiply(sd_I_c, cor_chol)';
  vector[I] I_h_0 = std_dev_h * sd_I_h; //individual-level fitness residual
  vector[I] I_s_0 = std_dev_s * sd_I_s; //individual-level fitness residual
  vector[I] I_m_0 = std_dev_m * sd_I_m; //individual-level fitness residual
  
  vector[nyear] Y_c_0 = std_dev_Y_c_0 * sd_Y_c_0; //year effects on clutch size
  vector[nyear] Y_h_0 = std_dev_Y_h_0 * sd_Y_h_0; //year effects on clutch size
  vector[nyear] Y_s_0 = std_dev_Y_s_0 * sd_Y_s_0; //year effects on clutch size
  vector[nyear] Y_m_0 = std_dev_Y_m_0 * sd_Y_m_0; //year effects on clutch size
}
model {
  //new objects declared in the model block are for efficient specification of
  //the model likelihood and will not be saved with the (transformed) parameters
  
  //separate individual intercepts and slopes for clutch size
  vector[I] I_c_mu = col(I_c, 1); //RN intercept
  vector[I] I_c_disc = col(I_c, 2); //variability
  
  //initialize linear predictors with fixed effects for phenotype and fitness
  vector[N] mu_c = X_clutch * coef_c;
  vector[N] mu_h = X_hatch * coef_h;
  vector[N] mu_s = mu_s_0 + X_survive * coef_s;
  vector[N] mu_m = mu_m_0 + X_mass * coef_m;
  vector[N] res_c; //variability on latent scale
  
  //add random effects to each model (excl. selection)
  for (n in 1 : N) {
    mu_c[n] += Y_c_0[year_id[n]] + I_c_mu[id[n]] ;
    res_c[n] = int_c + I_c_disc[id[n]];
    mu_h[n] += Y_h_0[year_id[n]] + I_h_0[id[n]];
    mu_s[n] += Y_s_0[year_id[n]] + I_s_0[id[n]] + O_s_0[n] * sd_O_s;
    mu_m[n] += Y_m_0[year_id[n]] + I_m_0[id[n]];
  }
  
  //add selection effects to fitness model
  for (n in 1 : N) {
    mu_h[n] +=   (b_h[1] + db_h[1] * dateWSD[n]) * I_c_mu[id[n]] + (b_h[2] + db_h[2] * dateWSD[n]) * I_c_disc[id[n]]
               + (q_h[1] + dq_h[1]) * square(I_c_mu[id[n]])
               + (q_h[2] + dq_h[2]) * square(I_c_disc[id[n]])
               + (q_h[3] + dq_h[3]) * (I_c_mu[id[n]] * I_c_disc[id[n]]);
  }
  
  for (n in 1 : N) {
    mu_s[n] +=   (b_s[1] + db_s[1] * dateWSD[n]) * I_c_mu[id[n]] + (b_s[2] + db_s[2] * dateWSD[n]) * I_c_disc[id[n]]
               + (q_s[1] + dq_s[1]) * square(I_c_mu[id[n]])
               + (q_s[2] + dq_s[2]) * square(I_c_disc[id[n]])
               + (q_s[3] + dq_s[3]) * (I_c_mu[id[n]] * I_c_disc[id[n]]);
  }
  
  for (n in 1 : N) {
    mu_m[n] +=   (b_m[1] + db_m[1] * dateWSD[n]) * I_c_mu[id[n]] + (b_m[2] + db_m[2] * dateWSD[n]) * I_c_disc[id[n]]
               + (q_m[1] + dq_m[1]) * square(I_c_mu[id[n]])
               + (q_m[2] + dq_m[2]) * square(I_c_disc[id[n]])
               + (q_m[3] + dq_m[3]) * (I_c_mu[id[n]] * I_c_disc[id[n]]);
  }
  
  //model likelihoods
  for (n in 1 : N) {
    target += cumulative_logit_lpmf(clutch[n] | mu_c[n], exp(res_c[n]), mu_c_0);
    target += cumulative_logit_lpmf(hatch[n] | mu_h[n],  exp(int_h), mu_h_0);
    target += binomial_logit_lpmf(survive[n] | hatch[n], mu_s[n]);
    target += normal_lpdf(mass[n] | mu_m[n], sd_sigma_m);
  }
  
  //priors
  
  //fixed effects
  int_c ~ normal(0,1);
  int_h ~ normal(0,1);
  mu_c_0 ~ normal(0, 1);
  mu_h_0 ~ normal(0, 1);
  mu_s_0 ~ normal(0, 1);
  mu_m_0 ~ normal(0, 1);
  coef_c ~ normal(0, 1);
  coef_h ~ normal(0, 1);
  coef_s ~ normal(0, 1);
  coef_m ~ normal(0, 1);
  b_h ~ normal(0, 1);
  q_h ~ normal(0, 1);
  b_s ~ normal(0, 1);
  q_s ~ normal(0, 1);
  b_m ~ normal(0, 1);
  q_m ~ normal(0, 1);
  db_h ~ normal(0, 1);
  dq_h ~ normal(0, 1);
  db_s ~ normal(0, 1);
  dq_s ~ normal(0, 1);
  db_m ~ normal(0, 1);
  dq_m ~ normal(0, 1);
  
  //random effects
  
  //SDs
  sd_sigma_m ~ cauchy(0, 1);
  sd_O_s ~ cauchy(0, 1);
  sd_I_c ~ cauchy(0, 1);
  sd_I_h ~ cauchy(0, 1);
  sd_I_s ~ cauchy(0, 1);
  sd_I_m ~ cauchy(0, 1);
  sd_Y_c_0 ~ cauchy(0, 1);
  sd_Y_h_0 ~ cauchy(0, 1);
  sd_Y_s_0 ~ cauchy(0, 1);
  sd_Y_m_0 ~ cauchy(0, 1);
  
  //z-scores
  O_s_0 ~ std_normal();
  to_vector(std_dev_c) ~ std_normal();
  to_vector(std_dev_h) ~ std_normal();
  to_vector(std_dev_s) ~ std_normal();
  to_vector(std_dev_m) ~ std_normal();
  to_vector(std_dev_Y_c_0) ~ std_normal();
  to_vector(std_dev_Y_h_0) ~ std_normal();
  to_vector(std_dev_Y_s_0) ~ std_normal();
  to_vector(std_dev_Y_m_0) ~ std_normal();
  
  //correlations
  cor_chol ~ lkj_corr_cholesky(2);
}
generated quantities {
  //reaction norm VCV
  matrix[ntrait, ntrait] R_I_c = cor_chol * cor_chol'; //correlation matrix
  matrix[ntrait, ntrait] S_I_c = diag_matrix(sd_I_c); //sd matrix
  matrix[ntrait, ntrait] P_I_c = S_I_c * R_I_c * S_I_c; //covariance matrix
}
