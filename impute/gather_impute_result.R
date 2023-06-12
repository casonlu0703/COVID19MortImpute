##########################################################################
##          Diagnostics for imputation results and shape data           ##
##########################################################################

rm(list = ls())

library(COVIDYPLL)
library(data.table)
library(ggplot2)
library(rstan)
library(MASS)
library(tidyverse)
library(brms)
library(ggpubr)
library(parallel)
library(loo)

i <-  8 # primary models used in the manuscript are model 8 (M1), 16 (M2), and 18 (M3)

fit_hurdle <- readRDS(paste0("impute/results/fit_hurdle_agg", i,".RDS"))
ix <- c(1:nrow(covid19d_cty))
ll <- extract_log_lik(fit_hurdle, "log_lik")
rm(fit_hurdle)
ll <- ll[, ix[!ix %in% impute_sample$ix_miss]]

loo_1 <- loo(ll, cores = 4)
# print(loo_1)
saveRDS(loo_1, paste0("impute/results/loo", i,".RDS"))

rm(list = ls())

i <- 8

loo_out <- readRDS(paste0("impute/results/loo", i,".RDS"))
print(loo_out)

ix_outlier <- which(loo_out$diagnostics$pareto_k > 0.7)

ix <- c(1:nrow(covid19d_cty))
ix_miss <- impute_sample$ix_miss
ix_nonmiss <- ix[!ix %in% ix_miss]
ix_outlier_o <- ix_nonmiss[ix_outlier]

outlier_dt <- covid19d_cty[ix_outlier_o,
                           .(fips, county_name, state, quarter, age_group, urban_rural_code, pop_size, covid_19_deaths)]
table(outlier_dt$age_group)
table(outlier_dt$urban_rural_code)
table(outlier_dt$quarter)
table(outlier_dt$quarter, outlier_dt$urban_rural_code)

write.csv(outlier_dt, paste0("impute/results/outlier from loo ", i,".csv"), row.names = F)


rm(list = ls())

i <- 18  # there are models 1, 2, 3 (4 and 5 are deleted)

if (i %in% c(1:2, 6)) {
  warmup <- 500 + 1
}

if (i %in% c(3:5, 7:26)) {
  warmup <- 1000 + 1
}

fit_hurdle <- readRDS(paste0("impute/results/fit_hurdle_agg", i,".RDS"))


g_stan_dig <- stan_diag(fit_hurdle)
ggarrange(g_stan_dig[[1]], g_stan_dig[[2]], g_stan_dig[[3]], nrow = 3)
ggsave(paste0("impute/results/stan_diag", i,".png"),
       device = "png", height = 10, width = 10)

pars <- c("bQ", "shape", "b_hu")
sum_estimates <- round(summary(fit_hurdle, pars = pars, probs = c(0.025, 0.975))$summary, 3)
saveRDS(sum_estimates, paste0("impute/results/sum_estimates_hurdle_agg", i,".RDS"))

traceplot(fit_hurdle, pars = pars)
ggsave(paste0("impute/results/traceplot of fixed effects_hurdle_agg", i,".png"), device = "png",
       height = 12, width = 14)


sampler_params <- get_sampler_params(fit_hurdle, inc_warmup = FALSE)
sampler_params_chain1 <- sampler_params[[1]]
colnames(sampler_params_chain1)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))

r_hat <- stan_rhat(fit_hurdle)
max_R_hat <- max(r_hat$data, na.rm = T)
fit_summary <- summary(fit_hurdle, probs = c(0.5))$summary
N <- dim(fit_summary)[[1]]
iter <- dim(rstan::extract(fit_hurdle)[[1]])[[1]]
pct_neff_ratio <- sum(fit_summary[,5] / iter < 0.001, na.rm = T)

diag_stats <- data.frame(diag_stats = c("Max Rhat", "Pct Ratio of N Eff to Sample Size <0.001"),
                         value = c(max_R_hat, pct_neff_ratio))
write.csv(diag_stats, paste0("impute/results/diag_stats", i,".csv"), row.names = F)


#### Sample data

ix_ymis <- grep("Ymi", names(fit_hurdle@sim$samples[[1]]))
ymis_draws <- mclapply(ix_ymis, function(x) {
  y <- c(fit_hurdle@sim$samples[[1]][[x]][warmup:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[2]][[x]][warmup:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[3]][[x]][warmup:fit_hurdle@sim$iter])
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
             ymis_draws = ymis_draws),
        paste0("impute/bayes_impute_agg", i,".RDS"))


# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
data(covid19d_cty)
covid19d_cty[, y := copy(covid_19_deaths)]
covid19d_cty[, suppress := ifelse(is.na(covid_19_deaths), "suppressed", "unsuppressed")]
covid19d_cty[, y_new := copy(y)]

covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ymis_draws[x, ]
  covid19d_cty$y_new[is.na(covid19d_cty$covid_19_deaths)] <- tmp_y

  covid19d_cty[, y_cate := as.character(y_new)]
  covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, quarter)]
  out[, pct := N / sum(N), by = .(quarter)]
  out[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "suppressed", "unsuppressed")]
  out[order(y_cate, quarter)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(mean_pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, quarter)]
covid19deaths_sum[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "imputed", "unsuppressed")]

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
ggsave(paste0("impute/results/impute distribution_hurdle_agg", i,".tiff"),
       device = "tiff", height = 12, width = 12, compression = "lzw", type = "cairo")


# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
exp_grid <- data.table(expand.grid(y_cate = c(paste0(c(0:19)), "20+"),
                                   age_group = levels(covid19d_cty$age_group),
                                   quarter = c(1:4)))

covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ymis_draws[x, ]
  covid19d_cty$y_new[is.na(covid19d_cty$covid_19_deaths)] <- tmp_y

  covid19d_cty[, y_cate := as.character(y_new)]
  covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, age_group, quarter)]
  out[, pct := N / sum(N), by = .(age_group, quarter)]
  out[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "suppressed", "unsuppressed")]
  out[order(y_cate, age_group)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(mean_pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, age_group, quarter)]
covid19deaths_sum[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "imputed", "unsuppressed")]
covid19deaths_sum <- merge(exp_grid, covid19deaths_sum,
                           by = c("y_cate", "age_group", "quarter"),
                           all.x = T)
covid19deaths_sum[, `:=` (quarter = paste0("Q", quarter),
                          age_group = paste0("Age ", age_group))]
covid19deaths_sum[, age_group := gsub("-", "\u2013", age_group)]

g2 <- ggplot(data = covid19deaths_sum[y_cate != "0"]) +
  geom_bar(aes(x = y_cate, y = mean_pct, fill = suppress),
           stat = "identity", color = "gray30", size = 0.5) +
  geom_errorbar(aes(x = y_cate, ymin = lb, ymax = ub),
                width = 0, color = "gray30", size = 1) +
  facet_grid(age_group ~ quarter, scale = "free") +
  xlab("\nNumber of COVID-19 deaths") +
  ylab("Proportion of counties\n") +
  ggtitle("Distribution of county COVID-19 deaths (positive only)") +
  scale_fill_manual(values = c("deepskyblue", "gray")) +
  theme_bw() +
  theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 18, colour = "gray20"),
        strip.text.y = element_text(size = 18, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.height = unit(0.3, 'in'),
        legend.key.width = unit(0.3, 'in'),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom")

ggsave(paste0("impute/results/impute distribution_age_hurdle_agg", i,".tiff"),
       plot = g2, device = "tiff", height = 10, width = 22, dpi = 300,
       units="in", compression = "lzw+p", type = "cairo")



# Calculate 95% credible intervals for the proportion of counties with 1-9 COVID-19 deaths
covid19deaths_dist <- mclapply(c(1:nrow(ymis_draws)), function(x) {
  tmp_y <- ymis_draws[x, ]
  covid19d_cty$y_new[is.na(covid19d_cty$covid_19_deaths)] <- tmp_y

  covid19d_cty[, y_cate := as.character(y_new)]
  covid19d_cty[, y_cate := ifelse(y_new >= 20, "20+", y_cate)]
  covid19d_cty$y_cate <- factor(covid19d_cty$y_cate, levels = paste0(c(0:19, "20+")))

  out <- covid19d_cty[, list(N = .N), by = .(y_cate, urban_rural_code, age_group, quarter)]
  out[, pct := N / sum(N), by = .(urban_rural_code, age_group, quarter)]
  out[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "suppressed", "unsuppressed")]
  out[order(y_cate, urban_rural_code)]
}, mc.cores = 6)

covid19deaths_dist <- rbindlist(covid19deaths_dist)
covid19deaths_sum <- covid19deaths_dist[, list(mean_pct = mean(pct),
                                               lb = quantile(pct, prob = 0.025),
                                               ub = quantile(pct, prob = 0.975)),
                                        by = .(y_cate, urban_rural_code, age_group, quarter)]
covid19deaths_sum[, suppress := ifelse(y_cate %in% paste0(c(1:9)), "imputed", "unsuppressed")]
covid19deaths_sum[, `:=` (quarter = paste0("Q", quarter),
                          age_group = paste0("Age ", age_group))]
covid19deaths_sum[, age_group := gsub("-", "\u2013", age_group)]


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
    xlab("\nNumber of COVID-19 Deaths") +
    ylab("Proportion of Counties\n") +
    ggtitle(paste0(x, " (positive only): number of counties = ", tmp_n_cnty)) +
    scale_fill_manual(values = c("deepskyblue", "gray")) +
    theme_bw() +
    theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
          strip.text.x = element_text(size = 18, colour = "gray20"),
          strip.text.y = element_text(size = 18, colour = "gray20", angle = 0),
          strip.background = element_rect(colour = NA, fill = "white"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.position = "bottom")
  ggsave(paste0("impute/results/impute distribution_hurdle_agg", i,"(", x, ").png"),
         plot = g2, device = "png", height = 10, width = 22, dpi = 300)
}



# posterior predictive check
ix_ysim <- grep("y_sim", names(fit_hurdle@sim$samples[[1]]))
ysim_draws <- mclapply(ix_ysim, function(x) {
  y <- c(fit_hurdle@sim$samples[[1]][[x]][warmup:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[2]][[x]][warmup:fit_hurdle@sim$iter],
         fit_hurdle@sim$samples[[3]][[x]][warmup:fit_hurdle@sim$iter])
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

covid19d_cty[, suppress := ifelse(is.na(covid_19_deaths), "suppressed", "unsuppressed")]
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

covid19deaths_dist <- mclapply(c(1:nrow(ysim_draws)), function(x) {
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
                          age_group = paste0("Age ", age_group))]
covid19deaths_sum[, age_group := gsub("-", "\u2013", age_group)]
covid19deaths_sum[, label := scales::percent(pct, accuracy = 2)]
covid19deaths_sum[, model := ifelse(i == 8, "M1",
                                    ifelse( i == 16, "M2",
                                            ifelse(i == 18, "M3", "NA")))]
saveRDS(covid19deaths_sum, paste0("impute/results/sim_data_dt", i, ".RDS"))

g2 <- ggplot(data = covid19deaths_sum[y_cate != "0"]) +
  geom_bar(aes(x = y_cate, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  facet_grid(age_group~quarter, scale = "free") +
  xlab("\nNumber of COVID-19 deaths") +
  ylab("Proportion of counties\n") +
  ggtitle("\n(B) Distribution of county-level COVID-19 deaths\n(positive only; simulated vs data)") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme(plot.title = element_text(size = 28, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 20, colour = "gray20"),
        strip.text.y = element_text(size = 20, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom")
# ggsave(paste0("impute/results/distribution simulated and data_hurdle_agg", i,".png"),
#        device = "png", height = 10, width = 18)


g1 <- ggplot(data = covid19deaths_sum[y_cate == "0"]) +
  geom_bar(aes(x = type, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  geom_text(aes(x = type, y = 0.5, label = label,
                group = type),
            stat = "identity", position = "stack", vjust = 1.5, size = 6) +
  facet_grid(age_group~quarter, scale = "free") +
  ylab("Proportion of counties\n") +
  ggtitle("\n(A) Proportion of counties with zero COVID-19 death\n(simulated vs data)") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme(plot.title = element_text(size = 28, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 20, colour = "gray20"),
        strip.text.y = element_text(size = 20, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom")
# ggsave(paste0("impute/results/distribution simulated and data prop zeros_hurdle_agg", i,".png"),
#        device = "png", height = 8, width = 10)

g3 <- ggarrange(g1, g2, nrow = 2, common.legend = T, legend = "bottom", heights = c(1, 1.1))
saveRDS(g3, paste0("impute/results/sim_data_plot", i, ".RDS"))

# ggsave(paste0("impute/results/distribution simulated and data prop hurdle_agg", i,".png"),
#        plot = g3, device = "png", height = 16, width = 14)
# ggsave(paste0("impute/results/distribution simulated and data prop hurdle_agg", i,".tiff"),
#        plot = g3, device = "tiff", height = 16, width = 14, dpi = 300,
#        units="in", compression = "lzw+p", type = "cairo")




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

for (j in unique(covid19d_cty$urban_rural_code)) {
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
  save_file <- paste0("impute/results/distribution simulated and data hurdle_agg", i,"(", j, ").png")
  ggsave(save_file, device = "png", height = 10, width = 18)
}


