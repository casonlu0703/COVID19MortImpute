rm(list = ls())

library(COVIDYPLL)
library(data.table)
library(ggplot2)
library(rstan)
library(MASS)
library(tidyverse)
library(brms)
library(ggpubr)

# ### Using brms to get the stan code to base on
# covid19d_cty$state_num <- as.numeric(as.factor(covid19d_cty$state))
# covid19d_cty$fips_num <- as.numeric(as.factor(covid19d_cty$fips))
# covid19d_cty[, l_pop_size := log(pop_size + 1)]
# covid19d_cty <- covid19d_cty[!is.na(covid_19_deaths)]
#
# data_ls <- list(
#   y = covid19d_cty$covid_19_death,
#   q2 = 1 * (covid19d_cty$quarter == 2), # quarters
#   q3 = 1 * (covid19d_cty$quarter == 3), # quarters
#   q4 = 1 * (covid19d_cty$quarter == 4), # quarters
#   state = covid19d_cty$state,
#   county = covid19d_cty$fips,
#   l_pop_size = covid19d_cty$l_pop_size
# )
# b_model <- bf(y ~ 1 + q2 + q3 + q4 + l_pop_size + (1 | state/county),
#               hu ~ 1 + q2 + q3 + q4 + (1 | state/county),
#               decomp = "QR")
#
#
# fit the model
# time_begin <- Sys.time()
# m2 <- brm(data = data_ls,
#           family = hurdle_gamma,
#           b_model,
#           prior = c(prior(normal(0,10), class = "Intercept"),
#                     prior(normal(0,1), class = "b"),
#                     prior(gamma(0.01, 0.01), class = "shape")),
#           iter = 10, chains = 1, cores = 1,
#           seed = 20200518)
# print(Sys.time() - time_begin)
# saveRDS(m2, "m2.RDS")
#
# stancode(m2)
#
# stan_data <- standata(m2)


#### Create stan model
covid19d_cty$state_num <- as.numeric(as.factor(covid19d_cty$state))
covid19d_cty$fips_num <- as.numeric(as.factor(covid19d_cty$fips))

covid19d_cty$urban_rural_code <- factor(covid19d_cty$urban_rural_code,
                                    levels = c("Noncore", "Medium metro", "Small metro",
                                               "Large fringe metro", "Micropolitan",
                                               "Large central metro"))
covid19d_cty$quarter <- factor(covid19d_cty$quarter)
covid19d_cty[, l_pop_size := log(pop_size + 1)]

X <- model.matrix( ~ quarter * age_group + quarter * urban_rural_code + l_pop_size, data = covid19d_cty)
X_hu <- model.matrix( ~ quarter * age_group + urban_rural_code, data = covid19d_cty)


covid19d_cty[, y := copy(covid_19_deaths)]

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
  M_4 = 1
)

begin_time <- Sys.time()
fit_hurdle <- stan(
  file = "inst/impute/stan/impute_hurdle_qr.stan",  # Stan program
  data = data_ls,         # named list of data
  chains = 3,             # number of Markov chains
  warmup = 500,           # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 3,              # number of cores (could use one per chain)
  refresh = 10,
  pars = c("bQ", "shape", "b_hu", "Ymi", "y_sim"),
  seed = 20210519
)
print(Sys.time() - begin_time)


saveRDS(fit_hurdle, "inst/impute/results/fit_hurdle_qr.RDS")


fit_hurdle <- readRDS("inst/impute/results/fit_hurdle_qr.RDS")

pars <- c("bQ", "shape", "b_hu")
sum_estimates <- round(summary(fit_hurdle, pars = pars, probs = c(0.025, 0.975))$summary, 3)
saveRDS(sum_estimates, "inst/impute/results/sum_estimates_hurdle_qr.RDS")

traceplot(fit_hurdle, pars = pars)
ggsave("inst/impute/results/traceplot of fixed effects_hurdle_qr.png", device = "png",
       height = 12, width = 14)

library(parallel)

ix_ymis <- grep("Ymi", names(fit_hurdle@sim$samples[[1]]))
ymis_draws <- mclapply(ix_ymis, function(x) {
  y <- c(fit_hurdle@sim$samples[[1]][[x]][501:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[2]][[x]][501:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[3]][[x]][501:fit_hurdle@sim$iter])
  round(y)
}, mc.cores = 6)

set.seed(521)
x_samp <- sample(c(1:length(ymis_draws[[1]])), 1000, rep = F)

ymis_draws <- mclapply(ymis_draws, function(x) {
  x[x_samp]
}, mc.cores = 6)
ymis_draws <- do.call(cbind, ymis_draws)

ix_miss <- which(is.na(covid19d_cty$covid_19_deaths))

saveRDS(list(ix_miss = ix_miss,
             ymis_draws = ymis_draws), "inst/impute/bayes_impute_qr.RDS")

# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
data(covid19d_cty)
covid19d_cty[, y := copy(covid_19_deaths)]
covid19d_cty[, suppress := ifelse(is.na(covid_19_deaths), "suppressed", "non-suppressed")]
covid19d_cty[, y_new := copy(y)]

covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ymis_draws[x, ]
  covid19d_cty$y_new[is.na(covid19d_cty$covid_19_deaths)] <- tmp_y

  covid19d_cty[, y_cate := as.character(y_new)]
  covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, quarter)]
  out[, pct := N / sum(N), by = .(quarter)]
  out[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "suppressed", "non-suppressed")]
  out[order(y_cate, quarter)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(mean_pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, quarter)]
covid19deaths_sum[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "imputed", "non-suppressed")]

g1 <- ggplot(data = covid19deaths_sum) +
  geom_bar(aes(x = y_cate, y = mean_pct, fill = suppress),
           stat = "identity", color = "gray30") +
  facet_wrap(.~quarter, scale = "free_y") +
  xlab("number of COVID-19 deaths") +
  ylab("proportion of counties") +
  ggtitle("(A) Distribution of COVID-19 deaths by quarters") +
  scale_fill_manual(values = c("deepskyblue", "gray")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 12, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")

g2 <- ggplot(data = covid19deaths_sum[y_cate != "0"]) +
  geom_bar(aes(x = y_cate, y = mean_pct, fill = suppress),
           stat = "identity", color = "gray30") +
  geom_errorbar(aes(x = y_cate, ymin = lb, ymax = ub),
                width = 0, color = "gray30", size = 0.5) +
  facet_wrap(.~quarter, scale = "free_y") +
  xlab("number of COVID-19 deaths") +
  ylab("proportion of counties") +
  ggtitle("(B) Distribution of COVID-19 deaths by quarters (death counts >0)") +
  scale_fill_manual(values = c("deepskyblue", "gray")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 12, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")
ggarrange(g1, g2, nrow = 2, common.legend = T)
ggsave("inst/impute/results/impute distribution_hurdle_qr.tiff", device = "tiff", height = 12, width = 12)


# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ymis_draws[x, ]
  covid19d_cty$y_new[is.na(covid19d_cty$covid_19_deaths)] <- tmp_y

  covid19d_cty[, y_cate := as.character(y_new)]
  covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, age_group, quarter)]
  out[, pct := N / sum(N), by = .(age_group, quarter)]
  out[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "suppressed", "non-suppressed")]
  out[order(y_cate, age_group)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(mean_pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, age_group, quarter)]
covid19deaths_sum[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "imputed", "non-suppressed")]
covid19deaths_sum[, `:=` (quarter = paste0("Q", quarter),
                          age_group = paste0("Age", age_group))]


g2 <- ggplot(data = covid19deaths_sum[y_cate != "0"]) +
  geom_bar(aes(x = y_cate, y = mean_pct, fill = suppress),
           stat = "identity", color = "gray30") +
  geom_errorbar(aes(x = y_cate, ymin = lb, ymax = ub),
                width = 0, color = "gray30", size = 0.5) +
  facet_grid(age_group ~ quarter, scale = "free") +
  xlab("Number of COVID-19 Deaths") +
  ylab("Proportion of Counties") +
  ggtitle("Distribution of county COVID-19 deaths (positive only)") +
  scale_fill_manual(values = c("deepskyblue", "gray")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 12, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")
# ggarrange(g1, g2, nrow = 2, common.legend = T)
ggsave("inst/impute/results/impute distribution_age_hurdle_qr.png",
       plot = g2, device = "png", height = 10, width = 18)

# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ymis_draws[x, ]
  covid19d_cty$y_new[is.na(covid19d_cty$covid_19_deaths)] <- tmp_y

  covid19d_cty[, y_cate := as.character(y_new)]
  covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, urban_rural_code, age_group, quarter)]
  out[, pct := N / sum(N), by = .(urban_rural_code, age_group, quarter)]
  out[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "suppressed", "non-suppressed")]
  out[order(y_cate, urban_rural_code)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(mean_pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, urban_rural_code, age_group, quarter)]
covid19deaths_sum[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "imputed", "non-suppressed")]
covid19deaths_sum[, `:=` (quarter = paste0("Q", quarter),
                          age_group = paste0("Age", age_group))]


rucc_ls <- unique(covid19deaths_sum$urban_rural_code)

master_dt <- data.table(expand.grid(y_cate = unique(covid19deaths_sum$y_cate),
                                    urban_rural_code = rucc_ls,
                                    age_group = unique(covid19deaths_sum$age_group),
                                    quarter = unique(covid19deaths_sum$quarter)))
covid19deaths_sum <- merge(master_dt, covid19deaths_sum,
                           by = c("y_cate", "urban_rural_code", "age_group", "quarter"),
                           all.x = T)
n_counties <- unique(covid19d_cty[, .(fips, urban_rural_code)])
n_counties <- n_counties[, list(N = .N), by = .(urban_rural_code)]


for (x in rucc_ls) {
  tmp_n_cnty <- n_counties[urban_rural_code == x]$N

  g2 <- ggplot(data = covid19deaths_sum[y_cate != "0" &
                                          urban_rural_code == x]) +
    geom_bar(aes(x = y_cate, y = mean_pct, fill = suppress),
             stat = "identity", color = "gray30") +
    geom_errorbar(aes(x = y_cate, ymin = lb, ymax = ub),
                  width = 0, color = "gray30", size = 0.5) +
    facet_grid(age_group ~ quarter, scale = "free") +
    xlab("Number of COVID-19 Deaths") +
    ylab("Proportion of Counties") +
    ggtitle(paste0(x, " (positive only): number of counties = ", tmp_n_cnty)) +
    scale_fill_manual(values = c("deepskyblue", "gray")) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          strip.text.x = element_text(size = 12, colour = "gray20"),
          strip.text.y = element_text(size = 12, colour = "gray20"),
          strip.background = element_rect(colour = NA, fill = "white"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "bottom")

  ggsave(paste0("inst/impute/results/impute distribution_hurdle_qr (", x, ").png"),
         plot = g2, device = "png", height = 10, width = 18)
}



# posterior predictive check
ix_ysim <- grep("y_sim", names(fit_hurdle@sim$samples[[1]]))
ysim_draws <- mclapply(ix_ysim, function(x) {
  y <- c(fit_hurdle@sim$samples[[1]][[x]][501:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[2]][[x]][501:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[3]][[x]][501:fit_hurdle@sim$iter])
  round(y)
}, mc.cores = 6)

set.seed(521)
x_samp <- sample(c(1:length(ysim_draws[[1]])), 1000, rep = F)

ysim_draws <- mclapply(ysim_draws, function(x) {
  x[x_samp]
}, mc.cores = 6)
ysim_draws <- do.call(cbind, ysim_draws)

# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
master_dt <- data.table(expand.grid(y_cate = paste0(c(0:19, "20+")),
                                    quarter = c(1:4),
                                    age_group = unique(covid19d_cty$age_group)))

covid19d_cty[, suppress := ifelse(is.na(covid_19_deaths), "suppressed", "non-suppressed")]
covid19d_cty[, y_new := copy(y)]
covid19d_cty[, y_cate := as.character(y_new)]
covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

sum_dist <- covid19d_cty[, list(N = .N), by = .(y_cate, quarter, age_group)]
sum_dist[, pct := N / sum(N), by = .(age_group, quarter)]
sum_dist[order(y_cate, age_group, quarter)]
sum_dist <- sum_dist[!is.na(y_cate)]

master_dt <- merge(master_dt, sum_dist, by = c("y_cate", "quarter", "age_group"), all.x = T)
# master_dt$pct[is.na(master_dt$pct)] <- 0
master_dt[, `:=` (type = "data", N = NULL)]

covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ysim_draws[x, ]
  covid19d_cty[, y_sim := tmp_y]

  covid19d_cty[, y_cate := as.character(y_sim)]
  covid19d_cty[, y_cate := ifelse(y_sim >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, quarter, age_group)]
  out[, pct := N / sum(N), by = .(quarter, age_group)]
  out[order(y_cate, quarter, age_group)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, quarter, age_group)]
covid19deaths_sum[, type := "simulated"]
covid19deaths_sum <- rbindlist(list(covid19deaths_sum, master_dt), use.names = T, fill = T)
covid19deaths_sum[, type := factor(type, levels = c("data", "simulated"))]
covid19deaths_sum[, `:=` (quarter = paste0("Q", quarter),
                          age_group = paste0("Age", age_group))]

ggplot(data = covid19deaths_sum[y_cate != "0"]) +
  geom_bar(aes(x = y_cate, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  facet_grid(age_group~quarter, scale = "free") +
  xlab("Number of COVID-19 Deaths") +
  ylab("Proportion of Counties") +
  ggtitle("Distribution of county-level COVID-19 deaths (positive only; simulated vs data)") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        strip.text = element_text(size = 12, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")
ggsave(paste0("inst/impute/results/distribution simulated and data_hurdle_qr.png"),
       device = "png", height = 10, width = 18)


covid19deaths_sum[, label := round(pct, 2)]

ggplot(data = covid19deaths_sum[y_cate == "0"]) +
  geom_bar(aes(x = type, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  geom_text(aes(x = type, y = 0.5, label = label,
                group = type),
            stat = "identity", position = "stack", vjust = 1.5, size = 4) +
  facet_grid(age_group~quarter, scale = "free") +
  ylab("Proportion of Counties") +
  ggtitle("Proportion of counties with zero COVID-19 death (simulated vs data)") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        strip.text = element_text(size = 12, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = "none")
ggsave(paste0("inst/impute/results/distribution simulated and data prop zeros_hurdle_qr.png"),
       device = "png", height = 8, width = 10)




# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths by RUCC
master_dt <- data.table(expand.grid(urban_rural_code = unique(covid19d_cty$urban_rural_code),
                                    age_group = unique(covid19d_cty$age_group),
                                    y_cate = paste0(c(0:19, "20+")),
                                    quarter = c(1:4)))

sum_dist <- covid19d_cty[, list(N = .N), by = .(y_cate, urban_rural_code, age_group, quarter)]
sum_dist[, pct := N / sum(N), by = .(urban_rural_code, age_group, quarter)]
sum_dist[order(urban_rural_code, y_cate, age_group, quarter)]
sum_dist <- sum_dist[!is.na(y_cate)]

master_dt <- merge(master_dt, sum_dist, by = c("urban_rural_code", "age_group", "y_cate", "quarter"), all.x = T)
master_dt[, `:=` (type = "data", N = NULL)]

covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ysim_draws[x, ]
  covid19d_cty[, y_sim := tmp_y]

  covid19d_cty[, y_cate := as.character(y_sim)]
  covid19d_cty[, y_cate := ifelse(y_sim >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(urban_rural_code, age_group, y_cate, quarter)]
  out[, pct := N / sum(N), by = .(urban_rural_code, age_group, quarter)]
  out[order(urban_rural_code, age_group, y_cate, quarter)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(urban_rural_code, age_group, y_cate, quarter)]
covid19deaths_sum[, type := "simulated"]
covid19deaths_sum <- rbindlist(list(covid19deaths_sum, master_dt), use.names = T, fill = T)
covid19deaths_sum[, type := factor(type, levels = c("data", "simulated"))]
covid19deaths_sum[, `:=` (quarter = paste0("Q", quarter),
                          age_group = paste0("Age", age_group))]

for (i in unique(covid19d_cty$urban_rural_code)) {
  tmp_n_cnty <- n_counties[urban_rural_code == i]$N
  ggplot(data = covid19deaths_sum[urban_rural_code == i & y_cate != "0"]) +
    geom_bar(aes(x = y_cate, y = pct, fill = type),
             stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
    facet_grid(age_group ~ quarter, scale = "free") +
    xlab("Number of COVID-19 Deaths") +
    ylab("Proportion of Counties") +
    ggtitle(paste0(i, ": number of counties = ", tmp_n_cnty," (positive only; simulated vs data)")) +
    scale_fill_manual(values = c("gray", "deepskyblue")) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          strip.text.x = element_text(size = 12, colour = "gray20"),
          strip.background = element_rect(colour = NA, fill = "white"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "bottom")
  save_file <- paste0("inst/impute/results/distribution simulated and data hurdle_qr (", i, ").png")
  ggsave(save_file, device = "png", height = 10, width = 18)
}






