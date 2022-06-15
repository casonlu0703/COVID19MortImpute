rm(list = ls())

library(COVIDYPLL)
library(data.table)
library(ggplot2)
library(viridis)
library(ggpubr)
library(dplyr)
library(openxlsx)

### Plots: by state, urban_rural_code, SVI category, and urban_rural_code plus SVI categories
death_samp <- bind_samples(impute_model = "m1")
year_rle <- NA

usa_ypll_ls <- calculate_usa_ypll(dt = death_samp, year_rle = year_rle)

death_samp[, urban_rural_svi := ifelse(!is.na(svi_cate),
                                       paste0(urban_rural_code, ";\n", svi_cate, " SVI"), NA)]
death_samp[, svi_cate := ifelse(!is.na(svi_cate), paste0(svi_cate, " SVI"), NA)]

var_vec <- c("state", "urban_rural_code", "svi_cate", "urban_rural_svi")
axis_text_size <- c(8, 10, 12, 10)
axis_text_angle <- c(0, 0, 0, 90)
usa_text_angle <- c(90, 0, 0, 0)
width_vec <- c(8, 12, 12, 8)
height_vec <- c(10, 12, 12, 8)

for (i in c(1:length(var_vec))) {

  byvar <- var_vec[i]
  g_aa_ypll_rate <- plot_age_adjusted_ypll(dt = death_samp, byvar = byvar,
                                           usa_ypll_ls = usa_ypll_ls, panel_letter = NULL,
                                           usa_text_angle = usa_text_angle[i],
                                           axis.text.x.size = axis_text_size[i],
                                           axis.text.x.angle = axis_text_angle[i])
  g_prop_by_age <- plot_prop_ypll_by_age(dt = death_samp, byvar = byvar,
                                         n_age_breaks = 4, year_rle = NA,
                                         usa_ypll_ls = usa_ypll_ls,
                                         panel_letter = NULL,
                                         axis.text.y.size = axis_text_size[i])
  g_ratio <- plot_ratio_ypll_per_death(dt = death_samp, byvar = byvar, year_rle = NA,
                                      usa_ypll_ls = usa_ypll_ls,
                                      panel_letter = NULL,
                                      axis.text.y.size = axis_text_size[i])


  # ggarrange(g_aa_ypll_rate$g_out,
  #           ggarrange(g_prop_by_age$g_out, g_ratio$g_out, widths = c(1, 0.75)),
  #           nrow = 2, heights = c(0.5, 1))
  # ggsave(file = paste0('inst/age_group_and_ratio (', byvar, ').png'),
  #        device = "png", width = width_vec[i], height = height_vec[i])

  ggsave(plot = g_aa_ypll_rate$g_out, paste0('inst/g_aa_ypll_rate (', byvar, ').png'),
         device = "png", width = 14, height = 6)

  ggsave(plot = g_prop_by_age$g_out, paste0('inst/g_prop_by_age (', byvar, ').png'),
         device = "png", width = 8, height = 12)

  ggsave(plot = g_ratio$g_out, paste0('inst/g_ratio (', byvar, ').png'),
         device = "png", width = width_vec[i], height = height_vec[i])


  wb <- createWorkbook()
  addWorksheet(wb, "Age-Standardized YPLL Rate")
  plot_data <- g_aa_ypll_rate$plot_data
  keep_cols <- grep("_mean", colnames(plot_data), value = T)
  keep_cols <- c(byvar, keep_cols)
  plot_data <- plot_data[, ..keep_cols]
  old_names <- c("pop_size_mean", "covid_19_deaths_mean", "covid19_death_rate_mean", "covid19_death_rate_aa_mean",
                 "tot_ypll_mean", "ypll_rate_mean", "ypll_rate_aa_mean")
  setcolorder(plot_data, c(byvar, old_names))
  setnames(plot_data, old_names, c("Population Size", "COVID-19 Deaths", "COVID-19 Death Rate", "Age-Standardized COVID-19 Death Rate",
                                   "Total YPLL", "YPLL Rate", "Age-Standardized YPLL Rate"))
  writeData(wb, "Age-Standardized YPLL Rate", plot_data, startRow = 3, startCol = 1)

  addWorksheet(wb, "Age Proportion")
  plot_data <- g_prop_by_age$plot_data
  keep_cols <- c(byvar, "age_group_plot", "ypll", "total_ypll", "prop_ypll_plot")
  plot_data <- plot_data[, ..keep_cols]
  setnames(plot_data, c("age_group_plot", "ypll", "total_ypll", "prop_ypll_plot"),
           c("Age Group", "Total YPLL", "Total YPLL by Category", "Proportion of Total YPLL Attributable to Age Group"))
  writeData(wb, "Age Proportion", plot_data, startRow = 3, startCol = 1)

  addWorksheet(wb, "Ratio YPLL per Death")
  plot_data <- g_ratio$plot_data
  keep_cols <- c(byvar, "covid_19_deaths", "tot_ypll", "ypll_per_death", "usa_ypll_per_death", "ratio")
  plot_data <- plot_data[, ..keep_cols]
  setnames(plot_data, c("covid_19_deaths", "tot_ypll", "ypll_per_death", "usa_ypll_per_death", "ratio"),
           c("COVID-19 Deaths", "Total YPLL",
             "YPLL per COVID-19 Death", "USA YPLL per COVID-19 Death",
             "Ratio of YPLL per Death by County Attributes Relative to YPLL per Death in USA"))
  writeData(wb, "Ratio YPLL per Death", plot_data, startRow = 3, startCol = 1)
  saveWorkbook(wb, file = paste0("inst/plot_data (", byvar, ").xlsx"), overwrite = T)
}



### Creating national table
usa_ypll_by_age <- calculate_ypll(death_samp, age_adjusted_output = F)
usa_ypll <- calculate_ypll(death_samp, age_adjusted_output = T, year_rle = year_rle)
usa_ypll[, `:=` (age_group = "Total")]

usa_ypll_by_age <- rbindlist(list(usa_ypll_by_age, usa_ypll), use.names = T, fill = T)

usa_ypll_by_age[, `:=` (covid_death = format(round(covid_19_deaths_mean), big.mark = ","),
                        covid_death_ci = paste0("(", format(round(covid_19_deaths_lb), big.mark = ","),
                                                "\U2012", format(round(covid_19_deaths_ub), big.mark = ","), ")"),
                        covid_death_rate = format(round(covid19_death_rate_mean, 2), big.mark = ","),
                        covid_death_rate_ci = paste0("(", format(round(covid19_death_rate_lb, 2), big.mark = ","),
                                                     "\U2012", format(round(covid19_death_rate_ub, 2), big.mark = ","), ")"),
                        covid_death_rate_aa = format(round(covid19_death_rate_aa_mean, 2), big.mark = ","),
                        covid_death_rate_aa_ci = paste0("(", format(round(covid19_death_rate_aa_lb, 2), big.mark = ","),
                                                        "\U2012", format(round(covid19_death_rate_aa_ub, 2), big.mark = ","), ")"),
                        tot_ypll = format(round(tot_ypll_mean), big.mark = ","),
                        tot_ypll_ci = paste0("(", format(round(tot_ypll_lb), big.mark = ","),
                                             "\U2012", format(round(tot_ypll_ub), big.mark = ","), ")"),
                        ypll_rate = format(round(ypll_rate_mean), big.mark = ","),
                        ypll_rate_ci = paste0("(", format(round(ypll_rate_lb), big.mark = ","),
                                              "\U2012", format(round(ypll_rate_ub), big.mark = ","), ")"),
                        ypll_rate_aa = format(round(ypll_rate_aa_mean, 2), big.mark = ","),
                        ypll_rate_aa_ci = paste0("(", format(round(ypll_rate_aa_lb, 2), big.mark = ","),
                                                 "\U2012", format(round(ypll_rate_aa_ub, 2), big.mark = ","), ")")
)]

death_var <- c("covid_death", "covid_death_rate", "covid_death_rate_aa",
               "tot_ypll", "ypll_rate", "ypll_rate_aa")
ci_vars <- paste0(death_var, "_ci")
out_var <- c("age_group", death_var, ci_vars)

out_table <- usa_ypll_by_age[, ..out_var]

out_table[, (ci_vars) := lapply(.SD, function(x) gsub("[[:space:]]", "", x)), .SDcols = ci_vars]


out_table <- melt(out_table, measure = lapply(death_var, function(x) c(x, paste0(x, "_ci"))),
                  value.name = death_var)

out_table[, variable := ifelse(variable == 1, "mean", "CI")]
out_table[, variable := factor(variable, levels = c("mean", "CI"))]
out_table <- out_table[order(age_group, variable)]

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", out_table, startRow = 3, startCol = 2)
saveWorkbook(wb, file = paste0("inst/death related outcome by age.xlsx"), overwrite = T)


#### Comparison between Predicted Death Counts and Death Counts from USA Facts Data

rm(list = ls())

library(COVIDYPLL)
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

data(covid19d_usafacts)
covid19d_usafacts[, data_type := "data"]

county_deaths_ls <- mclapply(c("m1", "m2", "m3"), function(x) {
  death_samp <- bind_samples(impute_model = x)
  year_rle <- NA

  cor_dt <- death_samp[, list(covid_19_deaths = sum(covid_19_deaths)), by = c("simno", "fips")]
  cor_dt <- merge(cor_dt, covid19d_usafacts[, .(fips, usafacts_death)],
                  by = "fips", all.x = T)

  cor_dt <- cor_dt[, list(cor = cor(covid_19_deaths, usafacts_death)**2), by = c("simno")]
  cor_dt <- cor_dt[, list(mean = mean(cor),
                          lb = quantile(cor, prob = 0.025),
                          ub = quantile(cor, prob = 0.975))]
  cor_dt[, impute_model := x]

  tmp_dt <- calculate_ypll(dt = death_samp, byvar = "fips",
                           year_rle = year_rle,
                           age_adjusted_output = T,
                           export_data_by_simno = T)

  tmp_dt$agg_sum[, impute_model := x]
  tmp_dt$sim_dt[, impute_model := x]

  return(list(tmp_dt$agg_sum, cor_dt, tmp_dt$sim_dt))
}, mc.cores = 3)

corr <- rbindlist(lapply(county_deaths_ls, `[[`, 2))
write.csv(corr, "inst/corr_coef_predict_vs_usafacts.csv", row.names = F)

county_ypll <- rbindlist(lapply(county_deaths_ls, `[[`, 1))


county_attr <- unique(covid19d_cty[, .(fips, county_name, state, urban_rural_code)])

death_dt <- lapply(county_deaths_ls, `[[`, 3)
death_dt <- lapply(death_dt, function(x) {
  x <- merge(x, covid19d_usafacts[, .(fips, usafacts_death)], by = "fips", all.x = T)
  x <- merge(x, county_attr[, .(fips, urban_rural_code)], by = "fips", all.x = T)
  x[, pct_diff := ifelse(usafacts_death == 0, 0,
                        (covid_19_deaths - usafacts_death) / usafacts_death * 100)]
  x[, list(median_pct_diff = median(pct_diff),
           lb_pct_diff = quantile(pct_diff, prob = 0.25),
           ub_pct_diff = quantile(pct_diff, prob = 0.75),
           sd_pct_diff = sd(pct_diff)), by = c("urban_rural_code")]
})

mod_ls <- c("M1", "M2", "M3")
death_dt <- lapply(c(1:3), function(x) {
  death_dt[[x]][, impute_model := mod_ls[x]]
  death_dt[[x]][, urban_rural_code := factor(urban_rural_code,
                                             levels = c("Large central metro", "Large fringe metro",
                                                        "Medium metro", "Small metro",
                                                        "Micropolitan", "Noncore"))]
  death_dt[[x]][order(urban_rural_code)]
})

death_dt <- rbindlist(death_dt)

ggplot(death_dt) +
  geom_errorbar(aes(xmin = lb_pct_diff, xmax = ub_pct_diff, y = urban_rural_code),
                width = 0) +
  geom_point(aes(x = median_pct_diff, y = urban_rural_code),
             shape = 21, stroke = 1.2, size = 2, fill = "white") +
  facet_wrap(~ impute_model, scales = "fixed", ncol = 3) +
  labs(x = "\nPercent Difference of Imputed Results Relative to County COVID-19 Death Counts\nReported in USA Facts",
       y = "County FIPS Code\n") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(colour = NA, fill = "white"),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank())

ggsave("inst/comparison between predicted to USAFacts.png", width = 12, height = 6, device = "png")


colset <- c("median_pct_diff", "lb_pct_diff", "ub_pct_diff", "sd_pct_diff")
death_dt[, (colset) := lapply(.SD, function(x) paste0(format(round(x, 2), nsmall = 2), "%")),
         .SDcols = colset]
write.csv(death_dt, "inst/difference between USAFacts and imputation.csv", row.names = F)


## pull SVI data
svi <- unique(covid19d_cty[, .(fips, svi_num)])
ccvi <- data.table(read.xlsx("inst/extdata/surgo_ccvi/ccvi.xlsx", sheet = "ccvi-US-county"))
ccvi[, ccvi := round(ccvi * 100, 2)]
setnames(ccvi, "FIPS", "fips")
vi_dt <- merge(svi, ccvi[, .(fips, ccvi)], by = c("fips"), all.x = T)
print(cor(vi_dt$svi_num, vi_dt$ccvi, use = "complete.obs"))

county_ypll <- merge(county_ypll, vi_dt, by = "fips", all.x = T)
county_ypll <- merge(county_ypll, covid19d_usafacts, by = "fips", all.x = T)

county_ypll[, list(cor_svi_rate_aa = cor(svi_num, covid19_death_rate_aa_mean, use = "complete.obs"),
                   cor_svi_rate = cor(svi_num, covid19_death_rate_mean, use = "complete.obs"),
                   cor_svi_cnt = cor(svi_num, covid_19_deaths_mean, use = "complete.obs"),
                   cor_svi_data = cor(svi_num, usafacts_rate, use = "complete.obs"),
                   cor_ccvi_rate_aa = cor(ccvi, covid19_death_rate_aa_mean, use = "complete.obs"),
                   cor_ccvi_rate = cor(ccvi, covid19_death_rate_mean, use = "complete.obs"),
                   cor_ccvi_cnt = cor(ccvi, covid_19_deaths_mean, use = "complete.obs"),
                   cor_ccvi_data = cor(ccvi, usafacts_rate, use = "complete.obs")),
            by = c("impute_model")] %>%
  write.csv("inst/correlation_svi_ccvi_covid_deaths.csv", row.names = F)

county_ypll[, svi_ptile := cut(svi_num, c(-Inf, seq(10, 100, 10)),
                               label = c("(0, 10]", "(10, 20]", "(20, 30]", "(30, 40]", "(40, 50]",
                                         "(50, 60]", "(60, 70]", "(70, 80]", "(80, 90]", "(90, 100]"))]

plot_dt0 <- county_ypll[!is.na(svi_num)]
plot_dt0 <- plot_dt0[, .(fips, covid19_death_rate_aa_mean, svi_ptile, impute_model)]
plot_dt0[, impute_model := toupper(impute_model)]
setnames(plot_dt0, "covid19_death_rate_aa_mean", "rate")

plot_dt1 <- county_ypll[!is.na(svi_num)]
plot_dt1 <- plot_dt1[, .(fips, usafacts_rate, svi_ptile)]
setnames(plot_dt1, "usafacts_rate", "rate")
plot_dt1[, impute_model := "Crude death rate from USAFacts"]
plot_dt <- rbindlist(list(plot_dt0, plot_dt1), use.name = T, fill = T)
setnames(plot_dt, "impute_model", "data_type")
plot_dt[, data_type := ifelse(data_type == "M1", "M1: Age-standaradized death rate",
                              ifelse(data_type == "M2", "M2: Age-standaradized death rate",
                                     ifelse(data_type == "M3", "M3: Age-standaradized death rate", data_type)))]

ggplot(plot_dt) +
  geom_smooth(method = stats::loess, se = F,
              aes(x = svi_ptile, rate,
                  color=data_type, group = data_type),
              position = position_dodge(width = 0.9), size = 1.5) +
  geom_boxplot(aes(x = svi_ptile, y = rate, color = data_type),
               outlier.shape = NA, size = 0.5, position = position_dodge(width = 0.9),
               alpha = 0) +
  scale_color_viridis(option = "D", begin = 0.3, end = 0.8, discrete = T) +
  scale_y_continuous(breaks = seq(0, 500, 100), limits = c(0, 500)) +
  labs(x = "\nCounty Social Vulnerability Index",
       y = "COVID-19 death rate per 100,000 population\n",
       title = "Correlation Between County Social Vulnerability and COVID-19 Death Rate") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("inst/S6 Fig - correlation between SVI and death rate.tiff",
       width = 10, height = 5, device = "tiff",
       compression = "lzw", type = "cairo")



county_attr <- unique(covid19d_cty[, .(fips, county_name, state, urban_rural_code)])
county_attr <- merge(county_attr, vi_dt, by = "fips", all.x = T)
county_attr <- merge(county_attr, covid19d_usafacts, by = "fips", all.x = T)

county_ypll <- merge(county_ypll, county_attr[, .(fips, county_name, state, urban_rural_code)],
                     by = "fips", all.x = T)
county_ypll <- county_ypll[order(impute_model, -covid19_death_rate_aa_mean)]

top_n <- 20

top_usafacts <- county_attr[order(-usafacts_rate)][1:top_n,
                                                   .(fips, county_name, state,
                                                     urban_rural_code, svi_num, ccvi,
                                                     usafacts_rate)][,
                                                                     impute_model := "USA Facts Data"]
top_usafacts[, usafacts_rate := round(usafacts_rate, 2)]
setnames(top_usafacts, "usafacts_rate", "statistics")

county_ypll[, statistics := paste0(round(covid19_death_rate_aa_mean, 2), "\n[",
                                   round(covid19_death_rate_aa_lb, 2), ", ",
                                   round(covid19_death_rate_aa_ub, 2), "]")]

top_m1 <- county_ypll[impute_model == "m1"][
  order(-covid19_death_rate_aa_mean)][1:top_n,
                                      .(impute_model, fips, county_name, state, urban_rural_code,
                                        svi_num, ccvi, statistics)]
top_m2 <- county_ypll[impute_model == "m2"][
  order(-covid19_death_rate_aa_mean)][1:top_n,
                                      .(impute_model, fips, county_name, state, urban_rural_code,
                                        svi_num, ccvi, statistics)]
top_m3 <- county_ypll[impute_model == "m3"][
  order(-covid19_death_rate_aa_mean)][1:top_n,
                                      .(impute_model, fips, county_name, state, urban_rural_code,
                                        svi_num, ccvi, statistics)]

top_dt <- rbindlist(list(top_m1, top_m2, top_m3))
top_dt <- rbindlist(list(top_usafacts, top_dt), use.names = T, fill = T)
top_dt[, county_state := paste0(county_name, ", ", state, " (", fips, ")")]
top_dt <- top_dt[, .(county_state, urban_rural_code, svi_num, ccvi, statistics, impute_model)]
write.csv(top_dt, "inst/top_ranked_county.csv", row.names = F)



#### Create maps

rm(list = ls())

library(COVIDYPLL)
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


death_samp <- bind_samples(impute_model = "m1")
year_rle <- NA

county_ypll <- calculate_ypll(dt = death_samp, byvar = "fips",
                              year_rle = year_rle,
                              age_adjusted_output = T,
                              export_data_by_simno = F)

fips_data <- unique(covid19d_cty[, .(fips, county_name, state)])
county_ypll_dt <- merge(fips_data, county_ypll, by = "fips")

write.xlsx(county_ypll_dt, file = paste0("inst/county_ypll.xlsx"), row.names = F)

tmp_brk <- round(quantile(county_ypll$ypll_rate_aa_mean,
                          prob = c(0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95), na.rm = T))
county_ypll[, `Age-standardized YPLL rate\nper 100,000 people` :=
              cut(ypll_rate_aa_mean, breaks = c(-Inf, tmp_brk, Inf),
                  labels = paste0(paste0(c(20, 40, 50, 60, 70, 80, 90, 95, ">95"),
                                         "-th percentile"), ": ",
                                  c("0-1304", "1304-1942", "1942-2250", "2250-2582", "2582-3023",
                                    "3023-3660", "3660-4896", "4896-6335", ">=6335")))]

tmp_brk <- round(quantile(county_ypll$covid19_death_rate_aa_mean,
                          prob = c(0.25, 0.5, 0.75)))
county_ypll[, `Age-standardized\nCOVID-19 death rate\nper 100,000 people` :=
         cut(covid19_death_rate_aa_mean,
             breaks = c(-Inf, tmp_brk, Inf),
             labels = paste0(c("low", "medium low", "medium high", "high")))]

data(covid19d_usafacts)
covid19d_usafacts[, data_type := "data"]

county_ypll <- merge(county_ypll, covid19d_usafacts[, .(fips, usafacts_rate)],
                     by = c("fips"), all.x = T)
# county_ypll[, `:=` (usafacts_rate.x = NULL, usafacts_rate.y = NULL)]

tmp_brk <- round(quantile(county_ypll$usafacts_rate,
                          prob = c(0.25, 0.5, 0.75)))
county_ypll[, `COVID-19 death rate\nper 100,000 people` :=
              cut(usafacts_rate,
                  breaks = c(-Inf, tmp_brk, Inf),
                  labels = paste0(c("low", "medium low", "medium high", "high")))]


# system.time(county_ypll_quarter <- calculate_ypll(dt = death_samp, byvar = c("fips", "quarter"),
#                                                   year_rle = year_rle,
#                                                   age_adjusted_output = TRUE, export_data_by_simno = TRUE))
# system.time(county_ypll_age <- calculate_ypll(dt = death_samp, byvar = c("fips", "quarter"),
#                                               year_rle = year_rle, age_adjusted_output = FALSE,
#                                               export_data_by_simno = TRUE))


## Code from here: https://mathewkiang.com/2017/01/16/using-histogram-legend-choropleths/
## and here: https://rud.is/b/2014/11/16/moving-the-earth-well-alaska-hawaii-with-r/
## Lower 48 states ----
lower_48 <-  c("28", "37", "40", "51", "54", "22", "26", "25", "16", "12", "31", "53", "35",
               "46", "48", "06", "01", "13", "42", "29", "08", "49", "47", "56", "36", "20",
               "32", "17", "50", "30", "19", "45", "33", "04", "11", "34", "24", "23",
               "10", "44", "21", "39", "55", "41", "38", "05", "18", "27", "09") #, "02", "15") # "02" Alaska; "15" Hawaii


##  Import US states and counties from 2018 CB shapefile ----
##  Data from: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
##  We could use map_data() -- but want this to be generalizable to all shp.
allcounties <- readOGR(dsn = "inst/extdata/cb_2018_us_county_500k",
                       layer = "cb_2018_us_county_500k")
allstates   <- readOGR(dsn = "inst/extdata/cb_2018_us_state_500k",
                       layer = "cb_2018_us_state_500k")

allcounties <- spTransform(allcounties,
                           CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
allstates <- spTransform(allstates,
                         CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))


## A little munging and subsetting for maps ----
allcounties@data$fips <- as.numeric(as.character(allcounties@data$GEOID))
allcounties@data$stateid <- as.character(allcounties@data$STATEFP)
allstates@data$stateid <- as.character(allstates@data$STATEFP)

## Alaska and Hawaii
akcounties <- allcounties[allcounties$STATEFP == "02", ]
akcounties <- elide(akcounties, rotate=-45)
akcounties <- elide(akcounties, scale=max(apply(bbox(akcounties), 1, diff)) / 2.3)
akcounties <- elide(akcounties, shift=c(-2400000, -2500000))
proj4string(akcounties) <- proj4string(allcounties)

alaska <- allstates[allstates$STATEFP == "02", ]
alaska <- elide(alaska, rotate=-45)
alaska <- elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska <- elide(alaska, shift=c(-2400000, -2500000))
proj4string(alaska) <- proj4string(allstates)

hicounties <- allcounties[allcounties$STATEFP == "15", ]
hicounties <- elide(hicounties, rotate=-35)
hicounties <- elide(hicounties, shift=c(5300000, -1700000)) # c(5000000, -1500000)
proj4string(hicounties) <- proj4string(allcounties)

hawaii <- allstates[allstates$STATEFP == "15", ]
hawaii <- elide(hawaii, rotate=-35)
hawaii <- elide(hawaii, shift=c(5300000, -1700000)) # c(5000000, -1500000)
proj4string(hawaii) <- proj4string(allstates)

## Only use lower 48 states
subcounties <- subset(allcounties, allcounties@data$state %in% lower_48)
substates <- subset(allstates, allstates@data$state %in% lower_48)

uscounties <- rbind(subcounties, akcounties, hicounties)
usstates <- rbind(substates, alaska, hawaii)

## Fortify into dataframes
uscounties_df <- tidy(uscounties, region = "GEOID")
uscounties_df$id <- as.numeric(uscounties_df$id)
usstates_df <- tidy(usstates, region = "GEOID")


map_death_rate_aa <- ggplot(data = county_ypll) +
  geom_map(aes(map_id = fips, fill = `Age-standardized\nCOVID-19 death rate\nper 100,000 people`),
           map = uscounties_df, color = "gray80", size = 0.01, alpha = 1)  +
  expand_limits(x = uscounties_df$long, y = uscounties_df$lat) +
  scale_fill_brewer(palette  = 'BuPu') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_path(data = usstates_df, aes(long, lat, group = group),
            color = "gray30", size = .3, alpha = 1) +
  coord_equal() +
  theme_map() +
  theme(plot.title = element_text(size = 24, face = "bold"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(title = "(B) Imputation results with provisional data:\nage-standardized COVID-19 death rate in 2020\n")


map_death_rate <- ggplot(data = county_ypll) +
  geom_map(aes(map_id = fips, fill = `COVID-19 death rate\nper 100,000 people`),
           map = uscounties_df, color = "gray80", size = 0.01, alpha = 1)  +
  expand_limits(x = uscounties_df$long, y = uscounties_df$lat) +
  scale_fill_brewer(palette  = 'BuPu') +
  # scale_fill_viridis(discrete = TRUE, direction = -1, alpha = 1, palette  = 'plasma') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_path(data = usstates_df, aes(long, lat, group = group),
            color = "gray30", size = .3, alpha = 1) +
  coord_equal() +
  theme_map() +
  theme(plot.title = element_text(size = 24, face = "bold"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(title = "(A) USAFacts data:\ncrude COVID-19 death rate in 2020\n")

g_map <- ggarrange(map_death_rate, map_death_rate_aa, nrow = 1)
ggsave(paste0("inst/impute/results/maps comparison (USAFacts vs M1).tiff"),
       plot = g_map, device = "tiff", height = 7, width = 22, dpi = 300,
       units="in", compression = "lzw+p", type = "cairo")



