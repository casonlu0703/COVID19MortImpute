// generated with brms 2.15.0
functions {
  /* hurdle gamma log-PDF of a single response
   * Args:
   *   y: the response value
   *   alpha: shape parameter of the gamma distribution
   *   beta: rate parameter of the gamma distribution
   *   hu: hurdle probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_gamma_lpdf(real y, real alpha, real beta, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             gamma_lpdf(y | alpha, beta);
    }
  }
  /* hurdle gamma log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args:
   *   y: the response value
   *   alpha: shape parameter of the gamma distribution
   *   beta: rate parameter of the gamma distribution
   *   hu: linear predictor for the hurdle part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_gamma_logit_lpdf(real y, real alpha, real beta, real hu) {
    if (y == 0) {
      return bernoulli_logit_lpmf(1 | hu);
    } else {
      return bernoulli_logit_lpmf(0 | hu) +
             gamma_lpdf(y | alpha, beta);
    }
  }
  // hurdle gamma log-CCDF and log-CDF functions
  real hurdle_gamma_lccdf(real y, real alpha, real beta, real hu) {
    return bernoulli_lpmf(0 | hu) + gamma_lccdf(y | alpha, beta);
  }
  real hurdle_gamma_lcdf(real y, real alpha, real beta, real hu) {
    return log1m_exp(hurdle_gamma_lccdf(y | alpha, beta, hu));
  }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=0> Nmi;  // number of missings
  int<lower=1> Jmi[Nmi];  // positions of missings
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_hu;  // number of population-level effects
  matrix[N, K_hu] X_hu;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_hu_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_hu_1;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_hu = K_hu - 1;
  matrix[N, Kc_hu] Xc_hu;  // centered version of X_hu without an intercept
  vector[Kc_hu] means_X_hu;  // column means of X_hu before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_hu) {
    means_X_hu[i - 1] = mean(X_hu[, i]);
    Xc_hu[, i - 1] = X_hu[, i] - means_X_hu[i - 1];
  }
}
parameters {
  // Missing data
  vector<lower=0.99, upper=9.01>[Nmi] Ymi;  // estimated missings
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
  vector[Kc_hu] b_hu;  // population-level effects
  real Intercept_hu;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  vector[N_3] r_3_hu_1;  // actual group-level effects
  vector[N_4] r_4_hu_1;  // actual group-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  r_3_hu_1 = (sd_3[1] * (z_3[1]));
  r_4_hu_1 = (sd_4[1] * (z_4[1]));
}
model {
  // likelihood including constants
  // vector combining observed and missing responses
  vector[N] Yl = Y;
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  // initialize linear predictor term
  vector[N] hu = Intercept_hu + Xc_hu * b_hu;
  Yl[Jmi] = Ymi;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    hu[n] += r_3_hu_1[J_3[n]] * Z_3_hu_1[n] + r_4_hu_1[J_4[n]] * Z_4_hu_1[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = shape * exp(-(mu[n]));
  }
  for (n in 1:N) {
    target += hurdle_gamma_logit_lpdf(Yl[n] | shape, mu[n], hu[n]);
  }
  // priors including constants
  target += normal_lpdf(b | 0, 10);
  target += normal_lpdf(Intercept | 0, 10);
  target += gamma_lpdf(shape | 0.01, 0.01);
  target += logistic_lpdf(Intercept_hu | 0, 1);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  target += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_2[1]);
  target += student_t_lpdf(sd_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_3[1]);
  target += student_t_lpdf(sd_4 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_4[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_hu_Intercept = Intercept_hu - dot_product(means_X_hu, b_hu);
  vector[N] y_sim;
  vector[N] mu_cal = Intercept + Xc * b;
  vector[N] hu_cal = Intercept_hu + Xc_hu * b_hu;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu_cal[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    mu_cal[n] = shape * exp(-(mu_cal[n]));
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    hu_cal[n] += r_3_hu_1[J_3[n]] * Z_3_hu_1[n] + r_4_hu_1[J_4[n]] * Z_4_hu_1[n];
  }
  for (n in 1:N) {
    if (bernoulli_rng(inv_logit(hu_cal[n]))) {
      y_sim[n] = 0;
    } else {
      y_sim[n] = gamma_rng(shape, mu_cal[n]);
    }
  }
}
