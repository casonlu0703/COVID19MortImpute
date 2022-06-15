functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=0> Nmi;  // number of missings
  int<lower=1> Jmi[Nmi];  // positions of missings
  int<lower=1> K;  // number of population-level effects / fixed-effects coefficients (including Intercept)
  matrix[N, K] X;  // population-level design matrix (fixed-effects variables, including Intercept)
  // data for group-level effects at the state level
  int<lower=1> N_1;  // number of state levels
  int<lower=1> M_1;  // number of coefficients per state
  int<lower=1> J_1[N];  // state grouping indicator per observation
  // state-level predictor values (all 1s)
  vector[N] Z_1_1;
  // data for group-level effects at the county level nested under states
  int<lower=1> N_2;  // number of county levels
  int<lower=1> M_2;  // number of coefficients per county
  int<lower=1> J_2[N];  // state grouping indicator per observation
  // county-level predictor values (all 1s)
  vector[N] Z_2_1;
  //int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1; // removing the intercept
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector<lower=2, upper=10>[Nmi] Ymi;  // estimated missings
  vector[Kc] b;  // fixed effect coefficients
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
}
model {
  // likelihood including constants
  // vector combining observed and missing responses
  vector[N] Yl = Y;
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  Yl[Jmi] = Ymi;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = shape * exp(-(mu[n]));
  }
  target += gamma_lpdf(Yl | shape, mu);

  // priors including constants
  target += normal_lpdf(b | 0, 1);
  target += normal_lpdf(Intercept | 0, 10);
  target += gamma_lpdf(shape | 0.01, 0.01);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  target += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  vector[N] y_sim;
  vector[N] mu_cal = Intercept + Xc * b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu_cal[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
  }
  for (n in 1:N) {
    // apply the inverse link function
    mu_cal[n] = shape * exp(-(mu_cal[n]));
  }
  for(n in 1:N) {
    y_sim[n] = gamma_rng(shape, mu_cal[n]);
  }
}
