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
  real hurdle_gamma_rng(real alpha, real beta, real hu) {
    if (bernoulli_rng(inv_logit(hu))) {
      return 0;
    } else {
      return gamma_rng(alpha, beta);
    }
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=0> Nmi;   // number of all missings
  int<lower=0> Nmi1;  // number of missings age 18-49
  int<lower=0> Nmi2;  // number of missings age 50-85+
  int<lower=1> Jmi1[Nmi1];  // positions of missings age 18-49
  int<lower=1> Jmi2[Nmi2];  // positions of missings age 50-85+
  int<lower=1> Jmi_miss1[Nmi1];  // positions among all missings age 18-49
  int<lower=1> Jmi_miss2[Nmi2];  // positions among all missings age 50-85+
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
  // state covid-19 deaths
  int<lower=1> n_state_d;
  vector[n_state_d] state_d;
  vector[n_state_d] sd_state_d;
  int<lower=1> n_gp;
  int<lower=1> gp_ix[N];
  // national covid-19 deaths
  int<lower=1> n_nat_d;
  vector[n_nat_d] nat_d;
  vector[n_nat_d] sd_nat_d;
  int<lower=1> n_gp_nat;
  int<lower=1> gp_nat_ix[N];
  real miss_sd1;
  real miss_sd2;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  // matrices for QR decomposition
  matrix[N, Kc] XQ;
  matrix[Kc, Kc] XR;
  matrix[Kc, Kc] XR_inv;
  int Kc_hu = K_hu - 1;
  matrix[N, Kc_hu] Xc_hu;  // centered version of X_hu without an intercept
  vector[Kc_hu] means_X_hu;  // column means of X_hu before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  // compute and scale QR decomposition
  XQ = qr_thin_Q(Xc) * sqrt(N - 1);
  XR = qr_thin_R(Xc) / sqrt(N - 1);
  XR_inv = inverse(XR);
  for (i in 2:K_hu) {
    means_X_hu[i - 1] = mean(X_hu[, i]);
    Xc_hu[, i - 1] = X_hu[, i] - means_X_hu[i - 1];
  }
}
parameters {
  // Missing data
  vector<lower=0.6, upper=9.4>[Nmi1] Ymi1;  // estimated missings for age 18-49
  vector<lower=0.6, upper=9.4>[Nmi2] Ymi2;  // estimated missings for age 50-85+
  vector[Kc] bQ;  // population-level effects
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
  vector[Nmi] Ymi;
  // initialize linear predictor term
  vector[N] mu = Intercept + XQ * bQ;
  // initialize linear predictor term
  vector[N] hu = Intercept_hu + Xc_hu * b_hu;

  Ymi[Jmi_miss1] = Ymi1;
  Ymi[Jmi_miss2] = Ymi2;
  mu = Intercept + XQ * bQ;
  hu = Intercept_hu + Xc_hu * b_hu;

  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += (sd_1[1] * z_1[1][J_1[n]]) * Z_1_1[n] + (sd_2[1] * z_2[1][J_2[n]]) * Z_1_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    hu[n] += (sd_3[1] * z_3[1][J_3[n]]) * Z_3_hu_1[n] + (sd_4[1] * z_4[1][J_4[n]]) * Z_4_hu_1[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = shape * exp(-(mu[n]));
  }
}
model {
  // likelihood including constants
  // vector combining observed and missing responses
  vector[N] Yl = Y;
  vector[N] p_pred; // prediction
  vector[N] y_pred; // prediction
  vector[n_gp] y_agg = rep_vector(0, n_gp); // aggregate y_pred
  vector[n_gp_nat] y_agg_nat = rep_vector(0, n_gp_nat); // aggregate y_pred at national level (7 values)
  real y_temp;
  real y_temp_nat;
  int tmp_gp;
  int tmp_gp_nat;

  p_pred = inv_logit(hu);
  y_pred = (shape ./ mu);
  y_pred = y_pred .* p_pred;

  // Aggregate predicted deaths by age group to state level
  for (j in 1:N) {
    tmp_gp = gp_ix[j];
    y_temp = y_pred[j];
    y_agg[tmp_gp] += y_temp;
  }

  // Aggregate predicted deaths by age group to national level
  for (j in 1:N) {
    tmp_gp_nat = gp_nat_ix[j];
    y_temp_nat = y_pred[j];
    y_agg_nat[tmp_gp_nat] += y_temp_nat;
  }

  Yl[Jmi1] = Ymi1;
  Yl[Jmi2] = Ymi2;
  for (n in 1:N) {
    target += hurdle_gamma_logit_lpdf(Yl[n] | shape, mu[n], hu[n]);
  }
  // priors including constants
  target += normal_lpdf(Ymi1 | 1, miss_sd1);
  target += normal_lpdf(Ymi2 | 1, miss_sd2);
  target += normal_lpdf(bQ | 0, 10);
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
  target += normal_lpdf(y_agg | state_d, sd_state_d);
  target += normal_lpdf(y_agg_nat | nat_d, sd_nat_d);
}
generated quantities {
  vector[N] y_sim;
  vector[N] Yl = Y;
  vector[N] log_lik;

  Yl[Jmi1] = Ymi1;
  Yl[Jmi2] = Ymi2;
  for (n in 1:N) {
    y_sim[n] = hurdle_gamma_rng(shape, mu[n], hu[n]);
    log_lik[n] = hurdle_gamma_logit_lpdf(Yl[n] | shape, mu[n], hu[n]);
  }
}

