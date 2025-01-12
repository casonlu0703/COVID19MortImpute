#########################################################################################
##          Aggregate imputed death results and generate figures and reports           ##
#########################################################################################

rm(list = ls())

library(ggplot2)
library(data.table)
library(parallel)
library(openxlsx)
library(rgdal)
library(rgeos)
library(viridis)
library(grid)
library(maptools)
library(broom)
library(mapproj)
library(ggthemes)
library(raster)
library(loo)
library(ggpubr)

fs <- paste0("impute/R/", list.files("impute/R"))
for(i in fs) {
  source(i)
}

#### Compare LOO among models 8, 15-20, 25
# vec_mod <- c(8, 15:20, 25)
# for (i in vec_mod) {
#   assign(paste0("loo", i), readRDS(paste0("impute/results/loo", i,".RDS")))
# }
#
# loo_compare(loo8, loo16)
# loo_compare(loo8, loo18)

#### Get state and national level results
i <- 18

# sum_est <- readRDS(paste0("impute/results/sum_estimates_hurdle_agg", i,".RDS"))

temp <- readRDS(paste0("impute/bayes_impute_agg", i,".RDS"))

ix_miss <- temp$ix_miss
ymis <- temp$ymis_draws

set.seed(20210621)

mort_all <- mort2020_old[state == "US", .(age_group, covid_19_deaths)]
mort_usa <- mort_all[, list(covid_19_deaths = sum(covid_19_deaths))]

death_dt <- mclapply(c(1:nrow(ymis)), function(x) {
  covid19d_cty_old$covid_19_deaths[ix_miss] <- ymis[x, ]
  sum_dt <- covid19d_cty_old[, list(covid_19_deaths = sum(covid_19_deaths, na.rm = T)),
                     by = .(fips, age_group, pop_size)]
  sum_dt[, list(sum_covid19_deaths = sum(covid_19_deaths)), by = .(age_group)][, sim := x]
}, mc.cores = 6)

death_dt <- rbindlist(death_dt)
death_dt <- merge(death_dt, mort_all, by = "age_group", all.x = T)

death_dt <- death_dt[, list(sum_covid19_deaths = mean(sum_covid19_deaths),
                            lb_d = quantile(sum_covid19_deaths, 0.025),
                            ub_d = quantile(sum_covid19_deaths, 0.975),
                            pct_diff = mean(round((sum_covid19_deaths - covid_19_deaths)/covid_19_deaths* 100, 2)),
                            pct_lb = quantile(round((sum_covid19_deaths - covid_19_deaths)/covid_19_deaths* 100, 2), 0.025),
                            pct_ub = quantile(round((sum_covid19_deaths - covid_19_deaths)/covid_19_deaths* 100, 2), 0.975)),
                     by = .(age_group)]

death_usa_dt <- mclapply(c(1:nrow(ymis)), function(x) {
  covid19d_cty_old$covid_19_deaths[ix_miss] <- ymis[x, ]
  sum_dt <- covid19d_cty_old[, list(covid_19_deaths = sum(covid_19_deaths, na.rm = T)),
                         by = .(fips, age_group, pop_size)]
  sum_dt[, list(sum_covid19_deaths = sum(covid_19_deaths))][, sim := x]
}, mc.cores = 6)

death_usa_dt <- rbindlist(death_usa_dt)
death_usa_dt[, covid_19_deaths := mort_usa[["covid_19_deaths"]]]
death_usa_dt <- death_usa_dt[, list(sum_covid19_deaths = mean(sum_covid19_deaths),
                            lb_d = quantile(sum_covid19_deaths, 0.025),
                            ub_d = quantile(sum_covid19_deaths, 0.975),
                            pct_diff = mean(round((sum_covid19_deaths - covid_19_deaths)/covid_19_deaths* 100, 2)),
                            pct_lb = quantile(round((sum_covid19_deaths - covid_19_deaths)/covid_19_deaths* 100, 2), 0.025),
                            pct_ub = quantile(round((sum_covid19_deaths - covid_19_deaths)/covid_19_deaths* 100, 2), 0.975))]
death_usa_dt[, age_group := "Total"]


death_all <- rbindlist(list(death_dt, death_usa_dt), use.names = TRUE, fill = TRUE)

death_all[, `:=` (covid19death = format(round(sum_covid19_deaths), big.mark = ","),
                  covid19death_ci = paste0("[", format(round(lb_d), big.mark = ","), ", ",
                                    format(round(ub_d), big.mark = ","), "]"),
                  pct_diff = paste0(round(pct_diff, 2), "%"),
                  pct_diff_ci = paste0("[", paste0(round(pct_lb, 2), "%"), ", ",
                                  paste0(round(pct_ub, 2), "%"), "]"))]
death_all <- death_all[, .(age_group, covid19death, covid19death_ci, pct_diff, pct_diff_ci)]

death_all <- melt(death_all, measure = list(c("covid19death", "covid19death_ci"),
                                            c("pct_diff", "pct_diff_ci")),
                  value.name = c("covid19deaths", "pct_diff"))
death_all[, variable := ifelse(variable == 1, "posterior mean", "credible CI")]
death_all[, variable := factor(variable, levels = c("posterior mean", "credible CI"))]
death_all <- death_all[order(age_group, variable)]

## State level distribution
state_sum_dt_all <- mclapply(c(1:nrow(ymis)), function(x) {
  covid19d_cty_old$covid_19_deaths[ix_miss] <- ymis[x, ]
  sum_dt <- covid19d_cty_old[, list(covid_19_deaths = sum(covid_19_deaths, na.rm = T)),
                         by = .(state, age_group)]
  sum_dt[, simno := x]
  return(sum_dt)
}, mc.cores = 6)

state_sum_dt_all <- rbindlist(state_sum_dt_all)
state_sum_dt <- state_sum_dt_all[, list(m = mean(covid_19_deaths),
                                    lb_d = quantile(covid_19_deaths, 0.025),
                                    ub_d = quantile(covid_19_deaths, 0.975)),
                     by = .(state, age_group)]

state_sum_dt <- merge(state_sum_dt, mort2020[, .(state, age_group, covid_19_deaths)],
                      by = c("state", "age_group"), all.x = T)

cor(state_sum_dt$m[state_sum_dt$age_group %in% c("18-29", "30-39")],
    state_sum_dt$covid_19_deaths[state_sum_dt$age_group %in% c("18-29", "30-39")])
plot(state_sum_dt$m[state_sum_dt$age_group %in% c("18-29", "30-39")],
     state_sum_dt$covid_19_deaths[state_sum_dt$age_group %in% c("18-29", "30-39")])


setnames(state_sum_dt, c("m", "covid_19_deaths"), c("predicted", "data"))

state_sum_dt_long <- melt(state_sum_dt, id.vars = c("state", "age_group"),
                          variable.name = "type",
                          measure.vars = c("predicted", "data"))

state_sum_dt[, type := "predicted"]
state_sum_dt[, hit_target := ifelse(data >= lb_d & data <= ub_d, TRUE, FALSE)]

state_sum_dt_long <- merge(state_sum_dt_long,
                           state_sum_dt[, .(state, age_group, type, lb_d, ub_d)],
                           by = c("state", "age_group", "type"),
                           all.x = T)
state_sum_dt[, type := NULL]
state_sum_dt[, ratio := predicted / data]


state_sum_dt_long$type <- factor(state_sum_dt_long$type, levels = c("predicted", "data"))
state_sum_dt_long[, age_group := gsub("-", "\u2013", age_group)]

g2 <- ggplot(data = state_sum_dt_long) +
  geom_errorbar(aes(x = age_group, ymin = lb_d, ymax = ub_d),
                width = 0.3, color = "gray30", alpha = 0.7) +
  geom_point(aes(x = age_group, y = value, color = type, shape = type, alpha = type),
             size = 1.5, stroke = 1, fill = "white") +
  scale_color_manual(values = c("gray30", "firebrick")) +
  scale_shape_manual(values = c(18, 21)) +
  scale_alpha_manual(values = c(0.7, 1)) +
  facet_wrap(~ state, scales = "free_y", ncol = 8) +
  ylab("Number of COVID-19 deaths") +
  xlab("Age group") +
  labs(title = "Number of COVID-19 Deaths by Age and State/Locality (Predicted vs Data)") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 10, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom")

ggsave(paste0("impute/results/S5 Fig - state_covid19d_dist_agg", i, ".tiff"),
       plot = g2, device = "tiff", height = 7, width = 12,
       compression = "lzw", type = "cairo")



## State distribution
setnames(state_sum_dt_all, "covid_19_deaths", "predicted")
state_sum_dt_all <- merge(state_sum_dt_all, mort2020[, .(state, age_group, covid_19_deaths)],
                      by = c("state", "age_group"), all.x = T)
setnames(state_sum_dt_all, "covid_19_deaths", "data")

corr <- paste0(round(cor(state_sum_dt_all$data, state_sum_dt_all$predicted), 3) * 100, "%")

# g1 <- ggplot(data = state_sum_dt_all) +
g1 <- ggplot(data = state_sum_dt) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.8,
              color = "gray30", linetype = 1, size = 0.8) +
  geom_point(aes(x = data, y = predicted), shape = 21,
             size = 1.5, stroke = 1, fill = "white", color = "royalblue") +
  scale_x_continuous(label = scales::comma) +
  scale_y_continuous(label = scales::comma) +
  labs(title = "Correlation Between Predicted and Observed\nAnnual COVID-19 Deaths by Age, State or Locality",
       caption = "Note: The predicted COVID-19 deaths are the predicted mean from 1,000 posterior samples.") +
  xlab("\nReported COVID-19 deaths by state/locality") +
  ylab("Predicted COVID-19 deaths by state/locality\n") +
  geom_label(
    label = paste0("x = y"),
    x = 8500,
    y = 7500,
    label.size = 0.35
  ) +
  geom_label(
    label = paste0("\u03C1", " = ", corr),
    x = 1000,
    y = 8500,
    label.size = 0,
    size = 5
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        strip.text = element_text(size = 12, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom")
saveRDS(g1, paste0("impute/results/corr_state_predicted_data", i, ".RDS"))


## State level statistics accuracy
state_sum_dt_all[, sqdiff := (predicted - data)^2] # square difference
state_sum_dt_all[, N := .N, by = .(simno)]
state_sum_dt_all[, rmse := sqrt(sum(sqdiff) / N),
                 by = .(simno)]
accuracy_measure_state <- unique(state_sum_dt_all[, .(simno, rmse, N)])

print("Accuracy measure")
accuracy_measure_state[, mean(rmse)]
accuracy_measure_state[, sd(rmse)]
mean(state_sum_dt$hit_target)

added_dt <- data.table(age_group = c("", "mean of accuracy measure by state",
                                     "sd of accuracy measure by state",
                                     "pct hit state-level targets"),
                       covid19deaths = c(NA, round(accuracy_measure_state[, mean(rmse)], 2),
                                              round(accuracy_measure_state[, sd(rmse)], 2),
                                              round(mean(state_sum_dt$hit_target), 2)))

death_all <- rbindlist(list(death_all, added_dt), use.name = T, fill = T)
write.xlsx(death_all, paste0("impute/results/covid19deaths_impute_age_group", i, ".xlsx"), row.names = F)




## National level distribution
nat_sum_dt_all <- mclapply(c(1:nrow(ymis)), function(x) {
  covid19d_cty_old$covid_19_deaths[ix_miss] <- ymis[x, ]
  sum_dt <- covid19d_cty_old[, list(covid_19_deaths = sum(covid_19_deaths, na.rm = T)),
                         by = .(age_group)]
  return(sum_dt)
}, mc.cores = 6)

nat_sum_dt_all <- rbindlist(nat_sum_dt_all)
nat_sum_dt <- nat_sum_dt_all[, list(m = mean(covid_19_deaths),
                                        lb_d = quantile(covid_19_deaths, 0.025),
                                        ub_d = quantile(covid_19_deaths, 0.975)),
                                 by = .(age_group)]

nat_sum_dt <- merge(nat_sum_dt, mort2020[state == "US", .(age_group, covid_19_deaths)],
                      by = c("age_group"), all.x = T)

nat_sum_dt_all <- merge(nat_sum_dt_all, mort2020[state == "US", .(age_group, covid_19_deaths)],
                        by = c("age_group"), all.x = T)

nat_sum_dt_all[, small := covid_19_deaths.x <= covid_19_deaths.y]
nat_sum_dt_all[, list(p = mean(small)), by = .(age_group)]

## Calculate annual distribution of the Bayesian simulation data

death_yr <- mclapply(c(1:nrow(ymis)), function(x) {
  covid19d_cty_old$covid_19_deaths[ix_miss] <- ymis[x, ]
  sum_dt <- covid19d_cty_old[, list(covid_19_deaths = sum(covid_19_deaths, na.rm = T)),
                         by = .(fips)]
  sum_dt[, cate := ifelse(covid_19_deaths >= 10, "10+", covid_19_deaths)]
  out_dt <- sum_dt[, list(sum_N = .N), by = .(cate)]
  out_dt[, pct := round(sum_N / sum(sum_N), 4)]
  out_dt[, cate := factor(cate, levels = c(0:9, "10+"))]
  out_dt <- out_dt[order(cate)]
}, mc.cores = 6)

death_yr <- rbindlist(death_yr)

death_yr[, list(m = round(mean(pct), 4)), by = .(cate)]


#### Bind figures
dt8 <- readRDS(paste0("impute/results/sim_data_dt8.RDS"))
dt16 <- readRDS(paste0("impute/results/sim_data_dt16.RDS"))
dt18 <- readRDS(paste0("impute/results/sim_data_dt18.RDS"))

dt <- rbindlist(list(dt8, dt16, dt18), use.names = T)

g1_1 <- ggplot(data = dt[y_cate != "0" & model == "M1"]) +
  geom_bar(aes(x = y_cate, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  facet_grid(age_group~quarter, scale = "free") +
  xlab("Number of COVID-19 deaths") +
  ylab("Proportion of counties") +
  ggtitle("M1") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 14, colour = "gray20", angle = 0),
        strip.text.y = element_text(size = 12, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 4.5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")

g1_2 <- ggplot(data = dt[y_cate != "0" & model == "M2"]) +
  geom_bar(aes(x = y_cate, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  facet_grid(age_group~quarter, scale = "free") +
  xlab("Number of COVID-19 deaths") +
  ylab("Proportion of counties") +
  ggtitle("M2") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 14, colour = "gray20", angle = 0),
        strip.text.y = element_text(size = 12, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 4.5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")

g1_3 <- ggplot(data = dt[y_cate != "0" & model == "M3"]) +
  geom_bar(aes(x = y_cate, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  facet_grid(age_group~quarter, scale = "free") +
  xlab("Number of COVID-19 deaths") +
  ylab("Proportion of counties") +
  ggtitle("M3") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 14, colour = "gray20", angle = 0),
        strip.text.y = element_text(size = 12, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 4.5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")

g1 <- ggarrange(g1_1, g1_2, g1_3, nrow = 3, common.legend = T,
                heights = c(1, 1, 1), align = "v")

df <- data.frame(x = 0.5, y = 0.5,
                 label = "(B) Distribution of county-level COVID-19 deaths\n(positive only; simulated vs data)")
g_tmp <- ggplot(df) +
  geom_point(aes(x = x, y = y), pch = 21, size = 0) +
  geom_text(aes(x = x, y = y, label = label), size = 6, fontface = "bold") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))

g1 <- ggarrange(g_tmp, g1, nrow = 2, heights = c(0.05, 0.95))

g2_1 <- ggplot(data = dt[y_cate == "0" & model == "M1"]) +
  geom_bar(aes(x = type, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  geom_text(aes(x = type, y = 0.5, label = label,
                group = type),
            stat = "identity", position = "stack", vjust = 1.5, size = 4) +
  facet_grid(age_group~quarter, scale = "free") +
  ylab("Proportion of counties") +
  xlab("Type") +
  ggtitle("M1") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 13, colour = "gray20", angle = 0),
        strip.text.y = element_text(size = 11, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = "bottom")


g2_2 <- ggplot(data = dt[y_cate == "0" & model == "M2"]) +
  geom_bar(aes(x = type, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  geom_text(aes(x = type, y = 0.5, label = label,
                group = type),
            stat = "identity", position = "stack", vjust = 1.5, size = 4) +
  facet_grid(age_group~quarter, scale = "free") +
  ylab("Proportion of counties") +
  xlab("Type") +
  ggtitle("M2") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 13, colour = "gray20", angle = 0),
        strip.text.y = element_text(size = 11, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = "bottom")


g2_3 <- ggplot(data = dt[y_cate == "0" & model == "M3"]) +
  geom_bar(aes(x = type, y = pct, fill = type),
           stat = "identity", color = "gray30", position = "dodge2", size = 0.3) +
  geom_text(aes(x = type, y = 0.5, label = label,
                group = type),
            stat = "identity", position = "stack", vjust = 1.5, size = 4) +
  facet_grid(age_group~quarter, scale = "free") +
  ylab("Proportion of counties") +
  xlab("Type") +
  ggtitle("M3") +
  scale_fill_manual(values = c("gray", "deepskyblue")) +
  theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 13, colour = "gray20", angle = 0),
        strip.text.y = element_text(size = 11, colour = "gray20", angle = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = "bottom")

g2 <- ggarrange(g2_1, g2_2, g2_3, nrow = 3, common.legend = T,
                heights = c(1, 1, 1), align = "v")

df <- data.frame(x = 0.5, y = 0.5,
                 label = "(A) Proportion of counties with zero COVID-19 death\n(simulated vs data)")
g_tmp <- ggplot(df) +
  geom_point(aes(x = x, y = y), pch = 21, size = 0) +
  geom_text(aes(x = x, y = y, label = label), size = 6, fontface = "bold") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))

g2 <- ggarrange(g_tmp, g2, nrow = 2, heights = c(0.05, 0.95))

g_bind <- ggarrange(g2, g1, widths = c(1, 1.2), ncol = 2)

ggsave("impute/results/S3 Fig - sim_data_plot all", plot = g_bind,
       device = "tiff", width = 16, height = 14,
       compression = "lzw", type = "cairo")



g8 <- readRDS(paste0("impute/results/corr_state_predicted_data8.RDS"))
g16 <- readRDS(paste0("impute/results/corr_state_predicted_data16.RDS"))
g18 <- readRDS(paste0("impute/results/corr_state_predicted_data18.RDS"))

df <- data.frame(x = 0.5, y = 0.5, label = "M1")
g_tmp <- ggplot(df) +
  geom_point(aes(x = x, y = y), pch = 21, size = 0) +
  geom_text(aes(x = x, y = y, label = label), size = 8, fontface = "bold") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))


g8 <- ggarrange(g_tmp, g8, nrow = 2, heights = c(0.05, 0.95))

df <- data.frame(x = 0.5, y = 0.5, label = "M2")
g_tmp <- ggplot(df) +
  geom_point(aes(x = x, y = y), pch = 21, size = 0) +
  geom_text(aes(x = x, y = y, label = label), size = 8, fontface = "bold") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))


g16 <- ggarrange(g_tmp, g16, nrow = 2, heights = c(0.05, 0.95))

df <- data.frame(x = 0.5, y = 0.5, label = "M3")
g_tmp <- ggplot(df) +
  geom_point(aes(x = x, y = y), pch = 21, size = 0) +
  geom_text(aes(x = x, y = y, label = label), size = 8, fontface = "bold") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))


g18 <- ggarrange(g_tmp, g18, nrow = 2, heights = c(0.05, 0.95))


ggarrange(g8, g16, g18, nrow = 1)
ggsave("impute/results/S4 Fig - corr_state_predicted_data.tiff",
       device = "tiff", width = 16, height = 5,
       compression = "lzw", type = "cairo")


#### CDR and ASDR by US census division

rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(parallel)
library(rgdal)
library(rgeos)
library(viridis)
library(grid)
library(maptools)
library(broom)
library(mapproj)
library(ggthemes)
library(openxlsx)

fs <- paste0("impute/R/", list.files("impute/R"))
for(i in fs) {
  source(i)
}

death_samp <- bind_samples(impute_model = "m1")
year_rle <- NA

## Create US census region and division: https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf
state_nm_dt <- data.table(state = state.abb,
                          state_name = state.name)
nyc_ls <- list(state = c("NYC", "DC"), state_name = c("New York City", "District of Columbia"))
state_nm_dt <- rbindlist(list(state_nm_dt, nyc_ls))

state_nm_dt[, region := ifelse(state_name %in% c("Connecticut", "Maine", "Massachusetts", "New Hampshire",
                                                 "Rhode Island", "Vermont", "New Jersey", "New York", "New York City", "Pennsylvania"), "Northeast",
                               ifelse(state_name %in% c("Indiana", "Illinois", "Michigan", "Ohio",
                                                        "Wisconsin", "Iowa", "Kansas", "Minnesota", "Missouri",
                                                        "Nebraska", "North Dakota", "South Dakota"), "Midwest",
                                      ifelse(state_name %in% c("Delaware", "District of Columbia", "Florida",
                                                               "Georgia", "Maryland", "North Carolina", "South Carolina",
                                                               "Virginia", "West Virginia", "Alabama", "Kentucky",
                                                               "Mississippi", "Tennessee", "Arkansas", "Louisiana",
                                                               "Oklahoma", "Texas"),"South",
                                             ifelse(state_name %in% c("Arizona", "Colorado", "Idaho", "New Mexico",
                                                                      "Montana", "Utah", "Nevada", "Wyoming",
                                                                      "Alaska", "California", "Hawaii", "Oregon",
                                                                      "Washington"), "West", "Other"))))]
state_nm_dt[, region := factor(region, levels = c("Northeast", "Midwest", "South", "West"))]

state_nm_dt[, division := ifelse(state_name %in% c("Connecticut", "Maine", "Massachusetts", "New Hampshire",
                                                   "Rhode Island", "Vermont"), "Div1: New England",
                                 ifelse(state_name %in% c("New Jersey", "New York", "New York City", "Pennsylvania"), "Div2: Middle Atlantic",
                                        ifelse(state_name %in% c("Indiana", "Illinois", "Michigan", "Ohio", "Wisconsin"), "Div3: East North Central",
                                               ifelse(state_name %in% c("Iowa", "Kansas", "Minnesota", "Missouri",
                                                                        "Nebraska", "North Dakota", "South Dakota"), "Div4: West North Central",
                                                      ifelse(state_name %in% c("Delaware", "District of Columbia", "Florida",
                                                                               "Georgia", "Maryland", "North Carolina", "South Carolina",
                                                                               "Virginia", "West Virginia"), "Div5: South Atlantic",
                                                             ifelse(state_name %in% c("Alabama", "Kentucky", "Mississippi",
                                                                                      "Tennessee"), "Div6: East South Central",
                                                                    ifelse(state_name %in% c("Arkansas", "Louisiana",
                                                                                             "Oklahoma", "Texas"), "Div7: West South Central",
                                                                           ifelse(state_name %in% c("Arizona", "Colorado", "Idaho", "New Mexico",
                                                                                                    "Montana", "Utah", "Nevada", "Wyoming"), "Div8: Mountain",
                                                                                  ifelse(state_name %in% c("Alaska", "California", "Hawaii", "Oregon",
                                                                                                           "Washington"), "Div9: Pacific", "Other")))))))))]

state_nm_dt[, division := factor(division, levels = c("Div1: New England", "Div2: Middle Atlantic", "Div3: East North Central",
                                                      "Div4: West North Central", "Div5: South Atlantic", "Div6: East South Central",
                                                      "Div7: West South Central", "Div8: Mountain", "Div9: Pacific", "Other"))]

state_nm_dt[, state_name := NULL]

death_samp <- merge(death_samp, state_nm_dt, by = "state", all = T)


county_ypll <- calculate_ypll(dt = death_samp, byvar = c("fips", "region", "division"),
                              year_rle = year_rle,
                              age_adjusted_output = T,
                              export_data_by_simno = F)

fips_data <- unique(covid19d_cty[, .(fips, county_name, state)])

data(covid19d_usafacts)

county_ypll <- merge(county_ypll, covid19d_usafacts[, .(fips, usafacts_rate)],
                     by = c("fips"), all.x = T)

col <- c("usafacts_rate", "covid19_death_rate_mean", "covid19_death_rate_aa_mean")

cty_m <- county_ypll[, lapply(.SD, mean), by = c("region", "division"), .SDcols = col][order(division)]
cty_m[, (col) := lapply(.SD, round), .SDcols = col]

cty_lb <- county_ypll[, lapply(.SD, quantile, prob = 0.25), by = c("division"), .SDcols = col]
cty_lb[, (col) := lapply(.SD, round), .SDcols = col]
setnames(cty_lb, col, paste0(col, "_lb"))

cty_ub <- county_ypll[, lapply(.SD, quantile, prob = 0.75), by = c("division"), .SDcols = col]
cty_ub[, (col) := lapply(.SD, round), .SDcols = col]
setnames(cty_ub, col, paste0(col, "_ub"))

cty_ci <- merge(cty_lb, cty_ub, by= "division", all = T)
cty_ci[, `:=` (usafacts_ci = paste0("[", usafacts_rate_lb, ", ", usafacts_rate_ub, "]"),
               covid19_death_rate_ci = paste0("[", covid19_death_rate_mean_lb, ", ", covid19_death_rate_mean_ub, "]"),
               covid19_death_rate_aa_ci = paste0("[", covid19_death_rate_aa_mean_lb, ", ", covid19_death_rate_aa_mean_ub, "]") )]

cty_ci <- cty_ci[, .(division, usafacts_ci, covid19_death_rate_ci, covid19_death_rate_aa_ci)]

cty_m <- merge(cty_m, cty_ci, by = "division")
cty_m <- cty_m[order(-usafacts_rate)][, CDR_usafacts_rank := c(1:.N)]
cty_m <- cty_m[order(-covid19_death_rate_mean)][, CDR_impute_rank := c(1:.N)]
cty_m <- cty_m[order(-covid19_death_rate_aa_mean)][, ASDR_impute_rank := c(1:.N)]

cty_m[, `:=` (CDR_usafacts = paste0(usafacts_rate, "\n", usafacts_ci),
              CDR_impute = paste0(covid19_death_rate_mean, "\n", covid19_death_rate_ci),
              ASDR_impute = paste0(covid19_death_rate_aa_mean, "\n", covid19_death_rate_aa_ci))]

col_order <- c("region", "division", "CDR_usafacts", "CDR_usafacts_rank",
               "CDR_impute", "CDR_impute_rank", "ASDR_impute", "ASDR_impute_rank")
cty_m <- cty_m[, ..col_order]
cty_m <- cty_m[order(division)]

write.xlsx(cty_m, "impute/results/rate by region and division.xlsx")





