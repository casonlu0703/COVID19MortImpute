rm(list = ls())

library(COVIDYPLL)
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


temp <- readRDS(paste0("inst/impute/bayes_impute_agg8.RDS"))
ymis8 <- temp$ymis_draws

temp <- readRDS(paste0("inst/impute/bayes_impute_agg16.RDS"))
ymis16 <- temp$ymis_draws

temp <- readRDS(paste0("inst/impute/bayes_impute_agg18.RDS"))
ymis18 <- temp$ymis_draws


mean8 <- data.table(t(t(apply(ymis8, 2, mean))))
mean8[, model := "M1"]

mean16 <- data.table(t(t(apply(ymis16, 2, mean))))
mean16[, model := "M2"]

mean18 <- data.table(t(t(apply(ymis18, 2, mean))))
mean18[, model := "M3"]

mean_dt <- rbindlist(list(mean8, mean16, mean18))
setnames(mean_dt, c("value", "model"))

ggplot(mean_dt, aes(x = value, fill = model, color = model)) +
  geom_density(alpha = 0.5, size = 0.3) +
  scale_fill_manual(values = c("tomato", "yellowgreen", "deepskyblue")) +
  scale_color_manual(values = c("tomato", "yellowgreen", "deepskyblue")) +
  ggtitle("Distribution of mean estimates across all\n23,353 suppressed COVID-19 death counts") +
  theme_classic() +
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
ggsave(paste0("inst/impute/results/dist_missing_mean.png"),
       device = "png", height = 5, width = 7)




ix_miss <- temp$ix_miss

miss_dt <- covid19d_cty[ix_miss]

set.seed(2021)

samp_set <- data.table(age_group = rep(levels(miss_dt$age_group), each = 5),
                       index = c(sort(sample(miss_dt[age_group == "18-29", which = T], 5, replace = F)),
                         sort(sample(miss_dt[age_group == "30-39", which = T], 5, replace = F)),
                         sort(sample(miss_dt[age_group == "40-49", which = T], 5, replace = F)),
                         sort(sample(miss_dt[age_group == "50-64", which = T], 5, replace = F)),
                         sort(sample(miss_dt[age_group == "65-74", which = T], 5, replace = F)),
                         sort(sample(miss_dt[age_group == "75-84", which = T], 5, replace = F)),
                         sort(sample(miss_dt[age_group == "85+", which = T], 5, replace = F))))


ymis8 <- data.table(ymis8[, samp_set$index])
setnames(ymis8, paste0("age ", samp_set$age_group, " missing id ", samp_set$index))
ymis8[, simno := c(1:1000)]
dt8 <- melt(ymis8, id.vars = "simno")
dt8[, model := "M1"]

ymis16 <- data.table(ymis16[, samp_set$index])
setnames(ymis16, paste0("age ", samp_set$age_group, " missing id ", samp_set$index))
ymis16[, simno := c(1:1000)]
dt16 <- melt(ymis16, id.vars = "simno")
dt16[, model := "M2"]

ymis18 <- data.table(ymis18[, samp_set$index])
setnames(ymis18, paste0("age ", samp_set$age_group, " missing id ", samp_set$index))
ymis18[, simno := c(1:1000)]
dt18 <- melt(ymis18, id.vars = "simno")
dt18[, model := "M3"]

dt <- rbindlist(list(dt8, dt16, dt18))

ggplot(dt, aes(x = value, fill = model, color = model)) +
  geom_density(alpha = 0.5, size = 0.3) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("tomato", "yellowgreen", "deepskyblue")) +
  scale_color_manual(values = c("tomato", "yellowgreen", "deepskyblue")) +
  theme_classic()
# dev.off()

