rm(list = ls())

# library(COVIDYPLL)
library(data.table)
library(rstan)
library(MASS)

# setwd("./impute")
# Sys.setenv("resnum" = 1)
resnum <- as.numeric(Sys.getenv("resnum"))

options(mc.cores = parallel::detectCores())

print(getwd())

# covid19d_cty <- get(load("data/covid19d_cty_old.rda"))
# mort2020 <- get(load("data/mort2020_old.rda"))

covid19d_cty <- get(load("data/covid19d_cty.rda"))
mort2020 <- get(load("data/mort2020.rda"))

#### Create stan model
covid19d_cty$state_num <- as.numeric(as.factor(covid19d_cty$state))
covid19d_cty$fips_num <- as.numeric(as.factor(covid19d_cty$fips))
covid19d_cty[, row_ix := c(1:.N)]

state_num <- unique(covid19d_cty[, .(state, state_num)])

covid19d_cty$urban_rural_code <- factor(covid19d_cty$urban_rural_code,
                                    levels = c("Noncore", "Medium metro", "Small metro",
                                               "Large fringe metro", "Micropolitan",
                                               "Large central metro"))
covid19d_cty$quarter <- factor(covid19d_cty$quarter)
covid19d_cty[, l_pop_size := log(pop_size + 1)]
covid19d_cty[, age_num := as.numeric(age_group)]

X <- model.matrix( ~ quarter * age_group + quarter * urban_rural_code + l_pop_size, data = covid19d_cty)
X_hu <- model.matrix( ~ quarter * age_group + urban_rural_code, data = covid19d_cty)

covid19d_cty[, y := copy(covid_19_deaths)]

## Get state covid-19 deaths counts by age group
us_mort2020 <- mort2020[state == "US", .(age_group, covid_19_deaths)]
mort2020 <- mort2020[state != "US", .(state, age_group, covid_19_deaths)] #, sd_covid19d)]
mort2020 <- merge(mort2020, state_num, by = "state", all.x = T)
mort2020[, age_num := as.numeric(age_group)]
mort2020[, gp_ix := c(1:.N)]
covid19d_cty <- merge(covid19d_cty, mort2020[, .(state, age_group, gp_ix)],
                      by = c("state", "age_group"), all.x = T)
covid19d_cty <- covid19d_cty[order(row_ix)]

X <- model.matrix( ~ quarter * age_group + quarter * urban_rural_code + age_group * urban_rural_code + l_pop_size, data = covid19d_cty)
X_hu <- model.matrix( ~ quarter * age_group + quarter * urban_rural_code + urban_rural_code, data = covid19d_cty)
sd_vec <- ifelse(mort2020$covid_19_deaths < 5, 1, mort2020$covid_19_deaths * 0.2)
sd_vec_nat <- us_mort2020$covid_19_deaths * 0.2
iters <- 4000
warmup <- 1000


if (resnum == 26) {
  ## Poisson assumption
  sd_vec <- ifelse(mort2020$covid_19_deaths < 1, 1, sqrt(mort2020$covid_19_deaths))
  sd_vec_nat <- sqrt(us_mort2020$covid_19_deaths)
  iters <- 5000
  warmup <- 1000
}

data_ls <- list(
  N = nrow(covid19d_cty),
  Y = ifelse(is.na(covid19d_cty$y), Inf, covid19d_cty$y),
  Jmi = which(is.na(covid19d_cty$y)),
  Nmi = length(which(is.na(covid19d_cty$y))),
  K = ncol(X),
  X = X,
  Z_1_1 = rep(1, nrow(covid19d_cty)),
  Z_2_1 = rep(1, nrow(covid19d_cty)),
  K_hu = ncol(X_hu),
  X_hu = X_hu,
  Z_3_hu_1 = rep(1, nrow(covid19d_cty)),
  Z_4_hu_1 = rep(1, nrow(covid19d_cty)),
  J_1 = covid19d_cty$state_num,
  J_2 = covid19d_cty$fips_num,
  J_3 = covid19d_cty$state_num,
  J_4 = covid19d_cty$fips_num,
  N_1 = length(unique(covid19d_cty$state_num)),
  M_1 = 1,
  N_2 = length(unique(covid19d_cty$fips_num)),
  M_2 = 1,
  N_3 = length(unique(covid19d_cty$state_num)),
  M_3 = 1,
  N_4 = length(unique(covid19d_cty$fips_num)),
  M_4 = 1,
  # State-level COVID-19 deaths
  state_d = mort2020$covid_19_deaths,
  sd_state_d = sd_vec,
  n_state_d = nrow(mort2020),
  n_gp = max(mort2020$gp_ix),
  # group indices
  gp_ix = covid19d_cty$gp_ix
)

data_ls$nat_d <- us_mort2020$covid_19_deaths
data_ls$sd_nat_d <- sd_vec_nat
data_ls$n_nat_d <- nrow(us_mort2020)
data_ls$n_gp_nat <- nrow(us_mort2020)
data_ls$gp_nat_ix <- covid19d_cty$age_num

if (resnum %in% 8) { # M1 in the manuscript
  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 15) {
  data_ls$miss_sd <- 5.0

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 16) { # M2 in the manuscript
  data_ls$miss_sd <- 10.0

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 17) {
  data_ls$miss_sd1 <- 5.0
  data_ls$miss_sd2 <- 10.0
  data_ls$Nmi1 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("18-29", "30-39", "40-49")])
  data_ls$Nmi2 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("50-64", "65-74", "75-84", "85+")])
  data_ls$Jmi1 <- which(is.na(covid19d_cty$y) & covid19d_cty$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi2 <- which(is.na(covid19d_cty$y) & covid19d_cty$age_group %in% c("50-64", "65-74", "75-84", "85+"))
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_miss1 <- which(tmp_dt$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi_miss2 <- which(tmp_dt$age_group %in% c("50-64", "65-74", "75-84", "85+"))

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior_bifur.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 18) { # M3 in the manuscript
  data_ls$miss_sd1 <- 5.0
  data_ls$miss_sd2 <- 20.0
  data_ls$Nmi1 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("18-29", "30-39", "40-49")])
  data_ls$Nmi2 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("50-64", "65-74", "75-84", "85+")])
  data_ls$Jmi1 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("18-29", "30-39", "40-49")))
  data_ls$Jmi2 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("50-64", "65-74", "75-84", "85+")))
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_miss1 <- which(tmp_dt$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi_miss2 <- which(tmp_dt$age_group %in% c("50-64", "65-74", "75-84", "85+"))

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior_bifur.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 19) {
  data_ls$miss_sd1 <- 6
  data_ls$miss_sd2 <- 15.0
  data_ls$Nmi1 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("18-29", "30-39", "40-49")])
  data_ls$Nmi2 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("50-64", "65-74", "75-84", "85+")])
  data_ls$Jmi1 <- which(is.na(covid19d_cty$y) & covid19d_cty$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi2 <- which(is.na(covid19d_cty$y) & covid19d_cty$age_group %in% c("50-64", "65-74", "75-84", "85+"))
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_miss1 <- which(tmp_dt$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi_miss2 <- which(tmp_dt$age_group %in% c("50-64", "65-74", "75-84", "85+"))

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior_bifur.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 20) {
  data_ls$miss_sd1 <- 7
  data_ls$miss_sd2 <- 15.0
  data_ls$Nmi1 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("18-29", "30-39", "40-49")])
  data_ls$Nmi2 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("50-64", "65-74", "75-84", "85+")])
  data_ls$Jmi1 <- which(is.na(covid19d_cty$y) & covid19d_cty$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi2 <- which(is.na(covid19d_cty$y) & covid19d_cty$age_group %in% c("50-64", "65-74", "75-84", "85+"))
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_miss1 <- which(tmp_dt$age_group %in% c("18-29", "30-39", "40-49"))
  data_ls$Jmi_miss2 <- which(tmp_dt$age_group %in% c("50-64", "65-74", "75-84", "85+"))

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior_bifur.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 21) { # structure the sd of the prior to be
  data_ls$n_miss_sd <- max(covid19d_cty$age_num)
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_gp <- tmp_dt$age_num
  data_ls$exp_mean <- 0.14

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_structure_prior.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "v_sd_miss", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 22) { # structure the sd of the prior to be
  data_ls$n_miss_sd <- max(covid19d_cty$age_num)
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_gp <- tmp_dt$age_num
  data_ls$exp_mean <- 0.5

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_structure_prior.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "v_sd_miss", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 23) { # structure the sd of the prior to be
  tmp_lvl <- lapply(levels(covid19d_cty$urban_rural_code),
                    function(x) paste0(x, "; ", levels(covid19d_cty$age_group)))
  tmp_lvl <- unlist(tmp_lvl)
  tmp_cate <- factor(paste0(covid19d_cty$urban_rural_code,
                            "; ", covid19d_cty$age_group),
                     levels = tmp_lvl)
  covid19d_cty[, au_gp_num := as.numeric(tmp_cate)]

  data_ls$n_miss_sd <- max(covid19d_cty$au_gp_num)
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_gp <- tmp_dt$au_gp_num
  data_ls$exp_mean <- 0.2

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_structure_prior.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "v_sd_miss", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% 24) { # structure the sd of the prior to be
  tmp_lvl <- lapply(levels(covid19d_cty$urban_rural_code),
                    function(x) paste0(x, "; ", levels(covid19d_cty$age_group)))
  tmp_lvl <- unlist(tmp_lvl)
  tmp_cate <- factor(paste0(covid19d_cty$urban_rural_code,
                            "; ", covid19d_cty$age_group),
                     levels = tmp_lvl)
  covid19d_cty[, au_gp_num := as.numeric(tmp_cate)]

  data_ls$n_miss_sd <- max(covid19d_cty$au_gp_num)
  tmp_dt <- covid19d_cty[is.na(y)]
  data_ls$Jmi_gp <- tmp_dt$au_gp_num
  data_ls$exp_mean <- 0.5

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_structure_prior.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "v_sd_miss", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
if (resnum %in% c(25, 26)) {
  data_ls$miss_sd1 <- 10.0
  data_ls$miss_sd2 <- 5.0
  data_ls$miss_sd3 <- 5.0
  data_ls$miss_sd4 <- 10.0
  data_ls$miss_sd5 <- 20.0
  data_ls$miss_sd6 <- 20.0
  data_ls$miss_sd7 <- 20.0

  data_ls$Nmi1 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("18-29")])
  data_ls$Nmi2 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("30-39")])
  data_ls$Nmi3 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("40-49")])
  data_ls$Nmi4 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("50-64")])
  data_ls$Nmi5 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("65-74")])
  data_ls$Nmi6 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("75-84")])
  data_ls$Nmi7 <- nrow(covid19d_cty[is.na(y) & age_group %in% c("85+")])

  data_ls$Jmi1 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("18-29")))
  data_ls$Jmi2 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("30-39")))
  data_ls$Jmi3 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("40-49")))
  data_ls$Jmi4 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("50-64")))
  data_ls$Jmi5 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("65-74")))
  data_ls$Jmi6 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("75-84")))
  data_ls$Jmi7 <- which(is.na(covid19d_cty$y & covid19d_cty$age_group %in% c("85+")))

  tmp_dt <- covid19d_cty[is.na(y)]

  data_ls$Jmi_miss1 <- which(tmp_dt$age_group %in% c("18-29"))
  data_ls$Jmi_miss2 <- which(tmp_dt$age_group %in% c("30-39"))
  data_ls$Jmi_miss3 <- which(tmp_dt$age_group %in% c("40-49"))
  data_ls$Jmi_miss4 <- which(tmp_dt$age_group %in% c("50-64"))
  data_ls$Jmi_miss5 <- which(tmp_dt$age_group %in% c("65-74"))
  data_ls$Jmi_miss6 <- which(tmp_dt$age_group %in% c("75-84"))
  data_ls$Jmi_miss7 <- which(tmp_dt$age_group %in% c("85+"))

  begin_time <- Sys.time()
  fit_hurdle <- stan(
    file = "stan/impute_hurdle_agg_nat_prior_7age.stan",  # Stan program
    # file = "inst/impute/stan/impute_hurdle_agg_nat_prior_7age.stan",  # Stan program
    data = data_ls,         # named list of data
    chains = 3,             # number of Markov chains
    warmup = warmup,           # number of warmup iterations per chain
    iter = iters,            # total number of iterations per chain
    cores = 3,              # number of cores (could use one per chain)
    refresh = 10,
    pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim", "log_lik"),
    seed = 20210519
  )
  print(Sys.time() - begin_time)
}
saveRDS(fit_hurdle, paste0("results/fit_hurdle_agg", resnum, "_new.RDS"))

