rm(list = ls())

library(COVIDYPLL)
library(ggplot2)
library(data.table)
library(viridis)
library(ggpubr)

setwd("./impute")

dat_ls <- data(package = "COVIDYPLL")
dat_ls <- dat_ls$results[, "Item"]
data(list = dat_ls, package = "COVIDYPLL")

covid19d_cty_old[, quarter := paste0("Q", quarter)]

# Calculate suppress values in each quarter
covid19d_cty_old[, suppress_covid_19_deaths := ifelse(is.na(covid_19_deaths), "suppressed: 1\u20139 COVID-19 deaths",
                                             ifelse(covid_19_deaths == 0, "unsuppressed: 0 COVID-19 deaths", "unsuppressed: >9 COVID-19 deaths"))]
covid19d_cty_old[, suppress_covid_19_deaths := factor(suppress_covid_19_deaths,
                                                  levels = c("unsuppressed: 0 COVID-19 deaths", "unsuppressed: >9 COVID-19 deaths", "suppressed: 1\u20139 COVID-19 deaths"))]

sum_suppress <- covid19d_cty_old[, list(N = .N), by = .(quarter, suppress_covid_19_deaths)]
sum_suppress[, pct := round(N / sum(N) * 100, 1), by = .(quarter)]
sum_suppress[, label := ifelse(pct < 5, NA, paste0(format(pct, nsmall = 1), "%"))]
sum_suppress$label[sum_suppress$quarter == 1 & sum_suppress$suppress_covid_19_deaths == "unsuppressed: >9 COVID-19 deaths"] <- NA
sum_suppress <- sum_suppress[order(quarter, suppress_covid_19_deaths)]

g1 <- ggplot(data = sum_suppress, aes(x = quarter, y = pct)) +
  geom_bar(aes(fill = suppress_covid_19_deaths),
           color = "gray30", stat = "identity", position = "stack", size = 0.5,
           alpha = 0.8) +
  geom_text(aes(x = quarter, y = pct, label = label,
                group = suppress_covid_19_deaths),
            stat = "identity", position = "stack", vjust = 1.5, size = 5) +
  scale_fill_manual(values = c("lightskyblue1", "deepskyblue2", "deepskyblue4")) +
  ylab("% of data elements") +
  xlab("Quarter") +
  ggtitle("(A) Distribution of suppressed and unsuppressed COVID-19 deaths") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "top")

county_agg_spread <- covid19d_cty_old[, list(covid_19_deaths = sum(covid_19_deaths)), by = .(fips, quarter)]
county_agg_spread[, dum_covid_case := ifelse(is.na(covid_19_deaths) | covid_19_deaths > 0, 1, 0)]
sum_spread <- county_agg_spread[, list(N = sum(dum_covid_case),
                                       pct = round(mean(dum_covid_case) * 100, 1)), by = .(quarter)]
sum_spread[, label := paste0(format(N, big.mark = ","), "\n(", pct, "%)")]

g2 <- ggplot(data = sum_spread, aes(x = quarter, y = N)) +
  geom_bar(fill = "royalblue", color = "gray30", stat = "identity",
           size = 0.5, alpha = 0.8) +
  geom_text(aes(x = quarter, label = label),
            vjust = -0.2, stat = "identity", size = 5) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 3500)) +
  ylab("Number of counties") +
  xlab("Quarter") +
  ggtitle("(B) Number of counties with >0 COVID-19 deaths") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        strip.text.x = element_text(size = 14, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom")

ggarrange(g1, g2, nrow = 2, heights = c(2, 1.8))

ggsave("results/pre_impute/S1 Fig-dist_supp_non_supp.tiff", height = 10, width = 10,
       device = "tiff", compression = "lzw", type = "cairo")



sum_suppress <- covid19d_cty_old[, list(N = .N), by = .(urban_rural_code, quarter, suppress_covid_19_deaths)]
sum_suppress[, pct := round(N / sum(N) * 100, 1), by = .(urban_rural_code, quarter)]
sum_suppress[, label := ifelse(pct < 5, NA, paste0(format(pct, nsmall = 1), "%"))]
sum_suppress <- sum_suppress[order(quarter, suppress_covid_19_deaths)]
sum_suppress[, urban_rural_code := factor(urban_rural_code,
                                          levels = c("Large central metro", "Large fringe metro",
                                                     "Medium metro", "Small metro",
                                                     "Micropolitan", "Noncore"))]

g1 <- ggplot(data = sum_suppress, aes(x = quarter, y = pct)) +
  geom_bar(aes(fill = suppress_covid_19_deaths),
           color = "gray30", stat = "identity", position = "stack", size = 0.5, alpha = 0.7) +
  geom_text(aes(x = quarter, y = pct, label = label,
                group = suppress_covid_19_deaths),
            stat = "identity", position = "stack", vjust = 1.5) +
  scale_fill_manual(values = c("lightskyblue1", "deepskyblue2", "deepskyblue4")) +
  facet_wrap(.~urban_rural_code, scale = "free") +
  ylab("% of data elements") +
  xlab("Quarter") +
  ggtitle("(B) Urban-rural code: distribution of suppressed and unsuppressed COVID-19 deaths") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 16, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "top")


# Calculate suppress values in each quarter by age group
sum_suppress <- covid19d_cty_old[, list(N = .N), by = .(quarter, age_group, suppress_covid_19_deaths)]
sum_suppress[, pct := round(N / sum(N) * 100, 1), by = .(quarter, age_group)]
sum_suppress[, label := ifelse(pct < 5, NA, paste0(format(pct, nsmall = 1), "%"))]
sum_suppress <- sum_suppress[order(quarter, suppress_covid_19_deaths)]
sum_suppress[, age_group := gsub("-", "\u2013", age_group)]

g2 <- ggplot(data = sum_suppress, aes(x = quarter, y = pct)) +
  geom_bar(aes(fill = suppress_covid_19_deaths),
           color = "gray30", stat = "identity", position = "stack", size = 0.5, alpha = 0.7) +
  geom_text(aes(x = quarter, y = pct, label = label,
                group = suppress_covid_19_deaths),
            stat = "identity", position = "stack", vjust = 1.5) +
  scale_fill_manual(values = c("lightskyblue1", "deepskyblue2", "deepskyblue4")) +
  facet_wrap(.~age_group, ncol = 4) +
  ylab("% of data elements") +
  xlab("Quarter") +
  ggtitle("(A) Age group: distribution of suppressed and unsuppressed COVID-19 deaths") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 16, colour = "gray20"),
        strip.background = element_rect(colour = NA, fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "top")


ggarrange(g2, g1, common.legend = T, nrow = 2)

ggsave("results/pre_impute/S2 Fig-dist_supp_non_supp_combine.tiff",
       height = 12, width = 12, unit = "in",
       device = "tiff", compression = "lzw", type = "cairo")

