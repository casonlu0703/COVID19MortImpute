###################################################
##       R functions for data preparation        ##
###################################################

#' @title Get 2000 US standard population
#' @param url A string of the url linking to the data of the 2000 US standard population.
#' @import data.table
#' @export
dl_us_standard_population <- function(url = "https://seer.cancer.gov/stdpopulations/stdpop.singleagesthru99.txt") {
  dldf <- RCurl::getURL(url)
  dldf <- strsplit(dldf, "\n")[[1]]
  tmp_l <- nchar(dldf[1])

  # data dictionary from: https://seer.cancer.gov/stdpopulations/stdpopdic.html
  dldf <- lapply(dldf, function(x) {
    col1 <- substr(x, 1, 3)
    col2 <- substr(x, 4, 6)
    col3 <- substr(x, 7, tmp_l)
    c(col1, col2, col3)
  })
  dldf <- data.table(do.call(rbind, dldf))
  setnames(dldf, c("standard", "age", "std_pop"))
  dldf <- dldf[standard == "205"]
  tmp_cols <-c("age", "std_pop")
  dldf[, (tmp_cols) := lapply(.SD, as.numeric), .SDcols = tmp_cols]
  dldf[, standard := NULL]
  dldf[, std_pop_million := round(std_pop / sum(std_pop) * 1000000)]

  age_breaks <- set_age_breaks()
  dldf[, age_group := cut(age, breaks = age_breaks$breaks, labels = age_breaks$labels, right = F)]
  dldf <- dldf[!is.na(age_group)]

  tmp_cols <- c("std_pop_million")
  std_pop_wgt <- dldf[, lapply(.SD, sum), by = .(age_group), .SDcols = tmp_cols]
  std_pop_wgt[, std_pop_wgt := std_pop_million / sum(std_pop_million)]
  std_pop_wgt[, age_group := factor(age_group, levels = set_age_breaks()$labels)]
  return(std_pop_wgt)
}

#' @title County level population size
#' @import data.table
#' @export
get_county_pop_size <- function() {
  tab <- read.csv("extdata/ACSST5Y2019.S0101_2021-05-11T113736/ACSST5Y2019.S0101_data_with_overlays_2021-05-11T113240.csv")
  tmp_names <- tab[1, ]
  col_ix <- c(1, 2, grep("Estimate!!Total", tmp_names))
  tab <- tab[, col_ix]

  tmp_names <- tmp_names[col_ix]
  tmp_names <- trimws(unlist(lapply(tmp_names, function(x) gsub("[^\\d]+", " ", x, perl = T))))
  tmp_names <- gsub(" ", "to", tmp_names)
  tmp_names <- ifelse(tmp_names != "", paste0("age", tmp_names), tmp_names)
  tmp_names[1:2] <- c("fips", "county_state")
  col_ix <- which(tmp_names != "")
  tmp_names <- unname(tmp_names[col_ix])

  tab <- tab[2:nrow(tab), col_ix]
  colnames(tab) <- tmp_names

  age_cols <- paste0("age", c("15to19", "15to17", "20to24", "25to29", "30to34", "35to39",
                              "40to44", "45to49", "50to54", "55to59", "60to64",
                              "65to69", "70to74", "75to79", "80to84", "85"))

  keep_cols <- c("fips", "county_state", age_cols)
  tab <- data.table(tab[, keep_cols])
  tab[, (age_cols) := lapply(.SD, as.numeric), .SDcols = age_cols]
  tab[, `:=` (age18to29 = (age15to19 - age15to17) + age20to24 + age25to29,
              age30to39 = age30to34 + age35to39,
              age40to49 = age40to44 + age45to49,
              age50to64 = age50to54 + age55to59 + age60to64,
              age65to74 = age65to69 + age70to74,
              age75to84 = age75to79 + age80to84,
              age85over = age85)]

  age_cols <- paste0("age", c("18to29", "30to39", "40to49", "50to64",
                              "65to74", "75to84", "85over"))
  keep_cols <- c("fips", "county_state", age_cols)

  tab <- unique(tab[, ..keep_cols])

  tab_long <- melt(tab, id.vars = c("fips", "county_state"), value.name = "pop_size")
  tab_long[, fips := as.numeric(gsub("0500000US", "", fips))]
  tmp_str <- do.call(rbind, strsplit(tab_long$county_state, ", "))
  tab_long[, `:=` (county_name = tmp_str[, 1],
                   state_name = tmp_str[, 2])]
  state_dt <- get_states()
  tab_long <- merge(tab_long, state_dt, by = c("state_name"), all.x = T)
  tab_long <- tab_long[!is.na(state)]

  # check whether there are 3142 counties
  n_cnty <- get_county_info()$n_cnty
  if (length(unique(tab_long$fips)) != n_cnty) {
    stop("The number of counties is not 3142")
  }

  tab_long[, age_group := gsub("to", "-", gsub("age", "", variable))]
  tab_long$age_group[tab_long$age_group == "85over"] <- "85+"
  tab_long[, age_group := factor(age_group, levels = set_age_breaks()$labels)]
  keep_cols <- c("fips", "county_name", "state", "age_group", "pop_size")
  tab_long <- tab_long[, ..keep_cols][order(fips, age_group)]
  tab_long$state[tab_long$fips %in% set_nyc_fips()] <- "NYC"
  return(tab_long)
}

#' @import data.table
#' @keywords internal
standardize_mort_data <- function(dt, state_level = FALSE) {
  dt <- data.table(dt)
  colnames(dt) <- tolower(colnames(dt))
  colnames(dt) <- gsub("\\.", "_", colnames(dt))

  tmp_start <- unlist(lapply(strsplit(dt$start_date, "/"), function(x) {
    paste0(x[3], "-", x[1], "-", x[2])
  }))
  tmp_end <- unlist(lapply(strsplit(dt$end_date, "/"), function(x) {
    paste0(x[3], "-", x[1], "-", x[2])
  }))

  dt[, `:=` (start_date = as.Date(tmp_start),
             end_date = as.Date(tmp_end))]

  dt$age_group[dt$age_group == "85 years and over"] <- "85+"
  dt[, age_group := gsub(" years", "", age_group)]
  if (isTRUE(state_level)) {
    dt[, age_group := factor(age_group, levels = c("All Ages", "0-17", set_age_breaks()$labels))]
  } else {
    dt[, age_group := factor(age_group, levels = set_age_breaks()$labels)]
  }
  dt <- dt[!is.na(age_group)]
  return(dt)
}

#' @title Provisional county-level COVID-19 deaths by quarter in 2020
#' @import data.table
#' @export
get_covid_death_cty <- function(version = "new") {
  if (version == "new") {
    tab <- read.csv("extdata/AH_Provisional_COVID-19_Deaths_by_Quarter__County_and_Age_for_2020 (2).csv")
  } else {
    tab <- read.csv("extdata/AH_Provisional_COVID-19_Deaths_by_Quarter__County_and_Age_for_2020.csv")
  }
  tab <- COVIDYPLL:::standardize_mort_data(tab)
  tab[, fips := copy(fips_code)]
  tab$state[tab$fips %in% set_nyc_fips()] <- "NYC"
  keep_cols <- c("fips", "county", "state", "urban_rural_code", "year", "quarter",
                 "start_date", "end_date", "age_group", "covid_19_deaths", "total_deaths")
  tab <- tab[, ..keep_cols][order(fips, quarter, age_group)]
  return(tab)
}

#' @title Provisional national and state-level COVID-19 and total deaths in 2020 and 2021
#' @import data.table
#' @export
get_mort_nation_state <- function(select_sex = "All Sexes", select_year = 2020, version = "new",
                                  verbose = TRUE) {
  if (version == "new") {
    tab <- read.csv("extdata/Provisional_COVID-19_Deaths_by_Sex_and_Age (2).csv")
  } else {
    tab <- read.csv("extdata/Provisional_COVID-19_Deaths_by_Sex_and_Age.csv")
  }
  tab <- standardize_mort_data(tab, state_level = T)
  tab <- tab[sex %in% select_sex & year %in% select_year & is.na(month)]
  state_dt <- get_states()
  state_dt <- rbindlist(list(state_dt, data.frame(state_name = "United States", state = "US")))
  tab <- merge(tab, state_dt, by.x = "state", by.y = "state_name", all.x = T)
  tab$state.y[tab$state == "New York City"] <- "NYC"
  tab[, state := NULL]
  setnames(tab, "state.y", "state")
  tab <- tab[!is.na(state)] # removing Puerto Rico

  if (isTRUE(verbose)) {
    print(paste0("There are ",
                 length(tab$covid_19_deaths[is.na(tab$covid_19_deaths) & !tab$age_group %in% c("0-17", "All Ages")]),
                 " rows of missing values. The total number of rows is ",
                 length(tab$covid_19_deaths[!tab$age_group %in% c("0-17", "All Ages")]),
                 " in the dataset.")) # 30 missing values among the dataset of 371 rows
  }

  ## Fill out NA values
  # assuming COVID-19 deaths in age 0-17 is ignorable (=0) if the value is NA
  tab_agg <- tab[age_group == "All Ages"]
  tab_agg <- tab_agg[, .(state, covid_19_deaths, total_deaths)]
  setnames(tab_agg, c("covid_19_deaths", "total_deaths"), c("tot_covid19d", "tot_d"))
  tab <- tab[age_group != "All Ages"]

  tab_agg[, `:=` (tot_covid19d = gsub(",", "", tot_covid19d),
                  tot_d = gsub(",", "", tot_d))]

  tab[, `:=` (covid_19_deaths = gsub(",", "", covid_19_deaths),
              total_deaths = gsub(",", "", total_deaths))]

  tab_agg[, `:=` (tot_covid19d = as.numeric(as.character(tot_covid19d)),
                  tot_d = as.numeric(as.character(tot_d)))]

  tab[, `:=` (covid_19_deaths = as.numeric(as.character(covid_19_deaths)),
              total_deaths = as.numeric(as.character(total_deaths)))]

  tab[, `:=` (n_covid19d = sum(covid_19_deaths, na.rm = T),
              n_death = sum(total_deaths, na.rm = T)),
      by = .(state)] #  get the sum of covid_19_deaths from the data

  tab <- merge(tab, tab_agg, by = c("state"), all.x = T)
  tab[, `:=` (remain_covid19d = tot_covid19d - n_covid19d)]

  tab[, `:=` (pr_covid19d = round(covid_19_deaths / n_covid19d, 3),
              pr_deaths = round(total_deaths / n_death, 3))]

  if (isTRUE(verbose)) {
    print(paste0("Corrrelation of age distribution between COVID-19 deaths and total deaths is ",
                 round(cor(tab$pr_covid19d[tab$age_group != "0-17"],
                           tab$pr_deaths[tab$age_group != "0-17"],
                           use = "complete.obs"), 3))) # the correlation is as high as 97.7%
  }

  # replace the missing values using the age distribution in total deaths
  # because there are no missing values among total deaths and
  # the correlation between covid19 deaths and total deaths is high
  tab[, tmp_group := ifelse(is.na(pr_covid19d), 1, 0)]
  tab[, renorm_pr_deaths := (tmp_group * pr_deaths) / sum((tmp_group * pr_deaths)),
      by = .(state)] # re-normalize the proportion of total deaths among the data rows that had missing COVID-19 deaths
  tab[, replace_missing_covid19d := round(renorm_pr_deaths * remain_covid19d)] # calculate the replaced values for the missing covid19 deaths using the remaining covid19 deaths times the renormalized proportion
  tab[, covid_19_deaths := ifelse(is.na(covid_19_deaths), replace_missing_covid19d, covid_19_deaths)] # replace missing COVID-19 deaths

  ## Drop All Ages and age group 0-17
  tab <- tab[age_group != "0-17"]
  tab$age_group <- droplevels(tab$age_group)

  keep_cols <- c("state", "year", "sex", "age_group", "covid_19_deaths", "total_deaths")
  tab <- tab[, ..keep_cols][order(state, age_group)]

  ## Get population estimates by state
  state_pop <- get_county_pop_size()
  state_pop$state[state_pop$fips %in% set_nyc_fips()] <- "NYC"
  state_pop <- state_pop[, list(pop_size = sum(pop_size)),
                         by = .(state, age_group)]

  tab <- merge(tab, state_pop, by = c("state", "age_group"), all.x = T)
  tab <- tab[order(state, age_group)]

  keep_cols <- c("state", "year", "sex", "age_group",
                 "covid_19_deaths", "total_deaths")
  tab <- tab[, ..keep_cols][order(state, age_group)]

  return(tab)
}

#' @title Postcensal population by age 2020
#' @import data.table
#' @export
get_postcensal_pop <- function() {
  tab <- read.csv("extdata/nc-est2019-alldata-r-file22.csv") # this is 2020 postcensal population
  tab <- data.table(tab)
  colnames(tab) <- tolower(colnames(tab))
  tab <- tab[month == 7 & (age > 17 & age < 900)] # using the population estimate in 7/1/2020
  tab <- tab[, .(age, tot_pop)]

  age_breaks <- set_age_breaks()
  tab[, age_group := cut(age, breaks = age_breaks$breaks, labels = age_breaks$labels, right = F)]
  tab <- tab[!is.na(age_group)]
  tab <- tab[, list(tot_pop = sum(tot_pop)), by = .(age_group)]
  return(tab)
}

#' @title Get provisional life expectancy estimates in 2017, 2018, and 2020
#' @import data.table
#' @export
calculate_provisional_le <- function() {
  ## calculate based on the technical note: https://www.cdc.gov/nchs/data/vsrr/VSRR10-508.pdf
  ## and here: https://www.statsdirect.com/help/survival_analysis/abridged_life_table.htm
  ## and here: https://www.cdc.gov/nchs/data/nvsr/nvsr61/nvsr61_03.pdf
  # mort2020 <- get_mort_nation_state()
  data(mort2020)
  age_breaks <- set_age_breaks()

  us_mort <- mort2020[state == "US", .(age_group, total_deaths)] # total deaths data
  us_pop <- get_postcensal_pop() # US population data

  us_mort <- merge(us_mort, us_pop, by = "age_group")
  us_mort[, nx := ifelse(is.infinite(diff(age_breaks$breaks)), 100 - 85, diff(age_breaks$breaks))]
  us_mort[, Mx := (total_deaths / tot_pop)]
  us_mort[, qx := (nx * Mx) / (1 + (1 - 0.5) * nx * Mx)] # assumed ax = 0.5
  us_mort[, qx := ifelse(qx > 1, 1, qx)]
  us_mort[, lx := 100000]
  us_mort[, dx := lx * qx]
  for (i in c(2:nrow(us_mort))) {
    us_mort$lx[i] <- us_mort$lx[i - 1] - us_mort$dx[i - 1]
    us_mort$dx[i] <- us_mort$lx[i] * us_mort$qx[i]
  }
  us_mort[, Lx := nx * (lx - dx) + 0.5 * nx * dx]
  us_mort$Lx[us_mort$age_group == "85+"] <- us_mort$lx[us_mort$age_group == "85+"] / us_mort$Mx[us_mort$age_group == "85+"]
  us_mort[, Tx := sum(Lx) - shift(cumsum(Lx), type = "lag")]
  us_mort$Tx[1] <- sum(us_mort$Lx)
  us_mort[, ex := Tx / lx]
  us_mort[, le2020 := copy(ex)]
  us_mort <- us_mort[, .(age_group, le2020)]

  ## Lift expectancy in 2018 https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/69-12/Table01.xlsx
  ## Lift expectancy in 2017 https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/68_07/Table01.xlsx
  wp <- c("2018" = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/69-12/Table01.xlsx",
          "2017" = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/68_07/Table01.xlsx")

  lftable <- lapply(c(1:length(wp)), function(i) {
    x <- wp[i]
    tmplft <- openxlsx::read.xlsx(x)
    colnames(tmplft) <- c("age", "qx", "lx", "dx", "Lx", "Tx", "ex")
    tmplft <- tmplft[3:(nrow(tmplft) - 1), ]
    tmplft <- data.table(tmplft)
    tmplft$age[tmplft$age == "100 and over"] <- "100–120"

    tmp_age <- strsplit(tmplft$age, "–")
    min_age <- as.numeric(unlist(lapply(tmp_age, `[[`, 1)))
    tmplft[, `:=` (min_age = min_age)]
    tmplft[, `:=` (age_group = cut(min_age, breaks = age_breaks$breaks,
                                    labels = age_breaks$labels, right = F),
                    age = NULL)]
    tmpcol <- c("qx", "lx", "dx", "Lx", "Tx", "ex")
    tmplft[, (tmpcol) := lapply(.SD, as.numeric), .SDcols = tmpcol]

    sublft <- tmplft[, list(max_lx = max(lx),
                                 max_Tx = max(Tx)),
                          by = .(age_group)]
    sublft <- sublft[, `:=` (min_lx = ifelse(is.na(shift(max_lx, type = "lead")), 0, shift(max_lx, type = "lead")),
                                     min_Tx = ifelse(is.na(shift(max_Tx, type = "lead")), 0, shift(max_Tx, type = "lead")))]

    sublft <- sublft[!is.na(age_group), ]
    sublft[, `:=` (dx = max_lx - min_lx,
                       Lx = max_Tx - min_Tx,
                       lx = max_lx)]
    sublft[, `:=` (qx = dx / lx,
                   Tx = max_Tx)]
    # sublft$Tx[1] <- sum(sublft$Lx)
    sublft[, paste0("le", names(x)) := Tx / lx]
    keep_cols <- c("age_group", grep("le", colnames(sublft), value = T))
    sublft[, ..keep_cols]
  })

  lftable <- Reduce(function(x, y) merge(x, y, by = "age_group", all.x = TRUE), lftable)
  lftable[, avg_le := (le2018 + le2017) / 2]

  us_mort <- merge(us_mort, lftable, by = "age_group", all.x = T)

  return(us_mort)
}

#' Get County-Level Percent Uninsured in 2019 (SAHIE estimates)
#' @import data.table
#' @export
get_uninsure_pop <- function() {
  tab <- read.csv("extdata/SAHIE_24JUN21_13_37_29_21.csv")
  keep_cols <- c("fips", "pct_uninsure")
  colnames(tab)[colnames(tab) %in% c("ID", "Uninsured...")] <- keep_cols
  tab <- data.table(tab)
  tab <- tab[, ..keep_cols]
  tab[, (keep_cols) := lapply(.SD, as.numeric), .SDcols = keep_cols]
  return(tab)
}

#' Get County-Level COVID-19 Deaths Data from USA Facts
#' @import data.table
#' @export
get_usafacts_covid_deaths <- function() {
  tab <- read.csv("extdata/covid_deaths_usafacts_20211118.csv")
  tab <- data.table(tab)
  setnames(tab, "countyFIPS", "fips")
  tab <- tab[fips != 0]
  keep_cols <- c("fips", grep("X20", colnames(tab), value = T))
  tab <- tab[, ..keep_cols]
  tab <- melt(tab, id.vars = "fips", variable.name = "date_0", value.name = "usafacts_death")
  tab[, date := as.Date(gsub("[.]", "-", gsub("X", "", date_0)))]
  tab[, date_0 := NULL]
  tab <- tab[order(fips, date)]
  tab <- tab[date == as.Date("2020-12-31")]

  pop_tab <- data.table(read.csv("extdata/covid_county_population_usafacts_20211118.csv"))
  setnames(pop_tab, "countyFIPS", "fips")

  tab <- merge(tab, pop_tab, by = "fips", all.x = T)
  tab[, usafacts_rate := usafacts_death / population * 100000]
  return(tab[, .(fips, usafacts_death, usafacts_rate)])
}
