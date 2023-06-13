###################################################
##                Create Data sets               ##
###################################################

rm(list = ls())

library(data.table)
library(parallel)

fs <- paste0("impute/R/", list.files("impute/R"))
for(i in fs) {
  source(i)
}

std_pop_wgt <- dl_us_standard_population()
le <- calculate_provisional_le() # life expectancy in 2017, 2018, and 2020
county_pop <- get_county_pop_size()
# covid19d_cty <- get_covid_death_cty()
covid19d_cty <- get_covid_death_cty(version = "old")
# mort2020 <- get_mort_nation_state()
mort2020 <- get_mort_nation_state(version = "old")
uninsure <- get_uninsure_pop()
covid19d_usafacts <- get_usafacts_covid_deaths()

## These imputation datasets are created after imputation is completed
impute_sample1 <- readRDS("inst/impute/bayes_impute_agg8.RDS") # aggregate results
impute_sample2 <- readRDS("inst/impute/bayes_impute_agg16.RDS") # aggregate results
impute_sample3 <- readRDS("inst/impute/bayes_impute_agg18.RDS") # aggregate results

# FIPS 2270 (Wade Hampton Census Area, AK) and and 46113 (Shannon County, SD)
# are not found in the census (county_pop) data
covid19d_cty[, `:=` (county = NULL, state = NULL)]
covid19d_cty <- merge(covid19d_cty, county_pop, by = c("fips", "age_group"), all.x = T)
# removing the two FIPS that had no population information
covid19d_cty <- covid19d_cty[!fips %in% c(2270, 46113)]
covid19d_cty[, row_ix := c(1:.N)]
covid19d_cty <- merge(covid19d_cty, uninsure, by = c("fips"), all.x = T)
covid19d_cty <- covid19d_cty[order(row_ix)]

# Policy data to get SVI categories
policy <- readRDS("extdata/policy_dt21.RDS")
keep_col <- c("fips", "quarter",
              "svi_overall", "svi_ses", "svi_household", "svi_minority", "svi_housing",
              "svi_overall_ter", "svi_ses_ter", "svi_household_ter", "svi_minority_ter", "svi_housing_ter")
policy <- policy[, ..keep_col]
setnames(policy, c(# "mask_acc1", "sah_acc1", "GB_acc1",
  "svi_overall", "svi_ses", "svi_household", "svi_minority", "svi_housing",
  "svi_overall_ter", "svi_ses_ter", "svi_household_ter", "svi_minority_ter", "svi_housing_ter"),
  c("svi_num", paste0("theme", c(1:4), "_num"),
    "svi_cate", paste0("theme", c(1:4), "_cate")))

covid19d_cty <- merge(covid19d_cty, unique(policy[, .(fips, svi_cate, svi_num)]), by = "fips", all.x = T)
covid19d_cty <- covid19d_cty[order(row_ix)]

policy[, svi_cate := NULL]

save(std_pop_wgt, file = "data/std_pop_wgt.rda")
save(le, file = "data/le.rda")
save(covid19d_cty, file = "data/covid19d_cty.rda")
# save(covid19d_cty_old, file = "data/covid19d_cty_old.rda")
save(mort2020, file = "data/mort2020.rda")
# save(mort2020_old, file = "data/mort2020_old.rda")
save(covid19d_usafacts, file = "data/covid19d_usafacts.rda")
save(impute_sample1, file = "data/impute_sample1.rda")
save(impute_sample2, file = "data/impute_sample2.rda")
save(impute_sample3, file = "data/impute_sample3.rda")
