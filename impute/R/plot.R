###################################################
##            R functions for figures            ##
###################################################

#' @title calculate national USA YPLL
#' @import data.table
#' @export
calculate_usa_ypll <- function(dt, year_rle = NA) {
  usa_ypll_by_age <- calculate_ypll(dt, age_adjusted_output = F)
  usa_ypll <- calculate_ypll(dt, age_adjusted_output = T, year_rle = year_rle)
  return(list(usa_ypll = usa_ypll, usa_ypll_by_age = usa_ypll_by_age))
}

#' @title bar charts
#' @import data.table
#' @export
plot_age_adjusted_ypll <- function(dt, byvar, year_rle = NA,
                      usa_ypll_ls = NULL,
                      panel_letter = NULL,
                      usa_text_angle = 90,
                      axis.text.x.size = 8, axis.text.x.angle = 30) {
  if (is.null(usa_ypll_ls)) {
    usa_ypll_ls = calculate_usa_ypll(dt, year_rle = year_rle)
  }
  usa_ypll <- usa_ypll_ls$usa_ypll
  usa_ypll[, paste0(byvar) := "USA"]

  sub_ypll <- calculate_ypll(dt, byvar = byvar, age_adjusted_output = T, year_rle = year_rle)

  all_ypll <- rbindlist(list(usa_ypll, sub_ypll), use.names = T)
  all_ypll <- all_ypll[order(-ypll_rate_aa_mean)]
  var_order <- unique(c("USA", as.character(all_ypll[, ..byvar][[byvar]])))
  all_ypll[, paste0(byvar) := factor(eval(parse(text = paste0(byvar))), levels = var_order)]
  all_ypll[, geo_level := ifelse(eval(parse(text = paste0(byvar))) == "USA", "national", byvar)]
  all_ypll[, text := ifelse(geo_level == "national", format(round(ypll_rate_aa_mean), big.mark = ","), NA)]
  all_ypll <- all_ypll[order(eval(parse(text = paste0(byvar))))]
  all_ypll <- all_ypll[!is.na(eval(parse(text = paste0(byvar))))]

  title_var <- ifelse(byvar == "state", "State and Locality",
                      ifelse(byvar == "svi_cate", "Overall Social Vulnerability",
                             ifelse(byvar == "urban_rural_code", "Urbanicity", "Urbanicity & Social Vulnerability")))
  if(!is.null(panel_letter)) panel_letter <- paste0(panel_letter, " ")
  title <- paste0(panel_letter, "Age-Standardized Years of Potential Life Lost Rate by ",
                    title_var, "\n(Per 100,000 Population)")

  hjust = -0.3
  vjust = 0.5
  if (usa_text_angle == 0) {
    hjust = 0.5
    vjust = -1
  }

  g_out <- ggplot(data = all_ypll) +
    geom_bar(aes(x = .data[[byvar]], y = ypll_rate_aa_mean, fill = geo_level),
             alpha = 0.9, color = "gray30", position = position_dodge(width = 0.9),
             stat = "identity", size = 0.3, vjust = vjust) +
    geom_text(aes(x = .data[[byvar]], y = ypll_rate_aa_mean, label = text, angle = usa_text_angle),
              size = 4, color = "black", hjust = hjust, vjust = vjust) +
    scale_y_continuous(labels = scales::comma) +
    facet_grid(.~geo_level, scales = "free", space="free") +
    ylab(paste0("Age-Standardized Years of Potential Life Lost Rate\n(Per 100,000 Population)")) +
    ggtitle(title) +
    scale_fill_manual(values = c("slateblue4", "deepskyblue3")) +
    theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          strip.text.x = element_blank(),
          panel.border = element_blank(),
          panel.spacing.x = unit(0.2, "line"),
          strip.background = element_rect(colour = NA, fill = "white"),
          legend.title = element_blank(),
          legend.text = element_blank(),
          axis.text.x = element_text(size = axis.text.x.size, angle = axis.text.x.angle, vjust = 0.5),
          axis.text.y = element_text(size = 10),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          legend.position = "none")
  return(list(g_out = g_out, plot_data = all_ypll))
}

#' @title Plot proportion of total YPLL attributable to each age group
#' @import data.table
#' @export
plot_prop_ypll_by_age <- function(dt, byvar,
                                  n_age_breaks = 7, # c(7, 4)
                                  year_rle = NA, usa_ypll_ls = NULL, panel_letter = NULL,
                                  axis.text.y.size = 8) {
  if (is.null(usa_ypll_ls)) {
    usa_ypll_ls <- calculate_usa_ypll(dt, year_rle = year_rle)
  }
  usa_ypll_by_age <- usa_ypll_ls$usa_ypll_by_age

  sub_ypll_by_age <- calculate_ypll(dt, byvar = byvar, age_adjusted_output = F)

  keep_col <- c("age_group", byvar, "covid_19_deaths_mean", "tot_ypll_mean")
  sub_prop <- sub_ypll_by_age[, ..keep_col]
  new_name <- c("covid_death", "ypll")
  setnames(sub_prop, c("covid_19_deaths_mean", "tot_ypll_mean"), new_name)

  sub_prop[, paste0("total_", new_name) := lapply(.SD, sum), by = c(byvar), .SDcol = new_name]
  sub_prop[, `:=` (prop_covid_death = round(covid_death / total_covid_death, 4),
                     prop_ypll = round(ypll / total_ypll, 4))]

  keep_col <- c("age_group", "covid_19_deaths_mean", "tot_ypll_mean")
  usa_prop <- usa_ypll_by_age[, ..keep_col]
  setnames(usa_prop, c("covid_19_deaths_mean", "tot_ypll_mean"), new_name)
  usa_prop[, paste0("total_", new_name) := lapply(.SD, sum), .SDcol = new_name]
  usa_prop[, `:=` (prop_covid_death = round(covid_death / total_covid_death, 4),
                   prop_ypll = round(ypll / total_ypll, 4))]
  usa_prop[, paste0(byvar) := "USA average"]

  age_key <- c("18-64", "18-64", "18-64", "18-64", "65-74", "75-84", "85+")
  names(age_key) <- set_age_breaks()$labels

  all_prop <- rbindlist(list(usa_prop, sub_prop), use.names = T)
  all_prop[, age_group2 := recode_factor(age_group, !!!age_key)]
  all_prop[, prop_ypll2 := sum(prop_ypll), by = c(byvar, "age_group2")]

  var_order <- as.character(all_prop[age_group == "18-29"][order(prop_ypll2)][[byvar]])
  var_order <- c("USA average", var_order[!var_order %in% "USA average"])

  all_prop[, paste0(byvar) := factor(eval(parse(text = paste0(byvar))), levels = var_order)]

  if (n_age_breaks == 4) {
    new_age_key <- c("18-39", "18-39", "40-64", "40-64", "65-84", "65-84", "85+")
    names(new_age_key) <- set_age_breaks()$labels
    all_prop[, age_group_plot := recode_factor(age_group, !!!new_age_key)]
  } else {
    all_prop[, age_group_plot := age_group]
  }

  all_prop <- all_prop[order(eval(parse(text = paste0(byvar))), age_group_plot)]
  all_prop[, `:=` (prop_ypll_plot = sum(prop_ypll),
                   ypll = sum(ypll)), by = c(byvar, "age_group_plot")]
  all_prop[, geo_level := ifelse(eval(parse(text = paste0(byvar))) == "USA average", "national", byvar)]
  age_order <- rev(levels(all_prop$age_group_plot))
  all_prop[, age_group_plot := factor(age_group_plot, levels = age_order)]
  all_prop[, text := ifelse(eval(parse(text = paste0(byvar))) == "USA average", prop_ypll_plot, NA)]
  all_prop <- all_prop[!is.na(eval(parse(text = paste0(byvar))))]

  title_var <- ifelse(byvar == "state", "State and Locality",
                      ifelse(byvar == "svi_cate", "Overall Social Vulnerability",
                             ifelse(byvar == "urban_rural_code", "Urbanicity", "Urbanicity & Social Vulnerability")))
  if(!is.null(panel_letter)) panel_letter <- paste0(panel_letter, " ")
  title <- paste0(panel_letter, title_var, ":\nProportion of Years of Potential Life Lost\nby Age Group")

  tmp_cols <- c("age_group_plot", "geo_level", "text", byvar, "ypll", "total_ypll", "prop_ypll_plot")
  uniq_all_prop <- unique(all_prop[, ..tmp_cols])
  age_color <- levels(uniq_all_prop$age_group_plot)[1:2]
  uniq_all_prop[, age_text_color := ifelse(age_group_plot %in% age_color & geo_level == "national", "black", "white")]
  setcolorder(uniq_all_prop, c(byvar, "geo_level", "age_group_plot"))

  g_out <- ggplot(data = uniq_all_prop) +
    geom_bar(aes(x = .data[[byvar]], y = prop_ypll_plot, fill = age_group_plot),
             color = "gray70", size = 0.1, position = "stack", stat="identity") +
    geom_text(aes(x = .data[[byvar]], y = prop_ypll_plot, label = text, color = age_text_color),
              position = position_stack(vjust = 0.5), size = 3) +
    facet_grid(geo_level~., scales = "free", space="free") +
    ggtitle(title) +
    ylab("Proportion of years of potential life lost") +
    coord_flip() +
    scale_color_manual(values = c("white", "black")) +
    scale_fill_viridis(begin = 0.1, end = 0.9,
                       discrete = T, option = "D", alpha = 0.9) +
    guides(colour = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          panel.border = element_blank(),
          panel.spacing.y = unit(0.2, "line"),
          strip.background = element_rect(colour = NA, fill = "white"),
          strip.text = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = axis.text.y.size),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.position = "right")
  return(list(g_out = g_out, plot_data = uniq_all_prop))
}

#' @title Plot ratio
#' @import data.table
#' @export
plot_ratio_ypll_per_death <- function(dt, byvar, year_rle = NA,
                                     usa_ypll_ls = NULL,
                                     panel_letter = NULL,
                                     axis.text.y.size = 8) {
  if (is.null(usa_ypll_ls)) {
    usa_ypll_ls <- calculate_usa_ypll(dt, year_rle = year_rle)
  }
  usa_ypll <- usa_ypll_ls$usa_ypll

  sub_ypll <- calculate_ypll(dt, byvar = byvar, age_adjusted_output = T, year_rle = year_rle)
  sel_cols <- c(byvar, "covid_19_deaths_mean", "tot_ypll_mean")
  ratio_ypll_death <- sub_ypll[, ..sel_cols]
  setnames(ratio_ypll_death, c("covid_19_deaths_mean", "tot_ypll_mean"),
           c("covid_19_deaths", "tot_ypll"))

  ratio_ypll_death[, `:=` (ypll_per_death = tot_ypll / covid_19_deaths,
                           usa_ypll_per_death = usa_ypll$tot_ypll_mean / usa_ypll$covid_19_deaths_mean)]

  ratio_ypll_death[, ratio := ypll_per_death / usa_ypll_per_death]
  ratio_ypll_death <- ratio_ypll_death[order(-ratio)]
  var_order <- as.character(ratio_ypll_death[, ..byvar][[byvar]])
  ratio_ypll_death[, paste0(byvar) := factor(eval(parse(text = paste0(byvar))), levels = rev(var_order))]
  ratio_ypll_death[, color := ifelse(ratio >= 1, "deepskyblue", "gray40")]

  ratio_ypll_death <- ratio_ypll_death[!is.na(eval(parse(text = paste0(byvar))))]

  title_var <- ifelse(byvar == "state", "State and Locality",
                      ifelse(byvar == "svi_cate", "Overall Social Vulnerability",
                             ifelse(byvar == "urban_rural_code", "Urbanicity", "Urbanicity & Social Vulnerability")))
  if(!is.null(panel_letter)) panel_letter <- paste0(panel_letter, " ")
  title <- paste0(panel_letter, "Ratio of YPLL per COVID-19 Death by ",  title_var, "\nRelative to YPLL per COVID-19 Death in USA")

  g_out <- ggplot(data = ratio_ypll_death) +
    geom_vline(xintercept = 1, color = "gray40", size = 0.3) +
    geom_point(aes(x = ratio, y = .data[[byvar]], color = color), stroke = 1) +
    geom_segment(aes(y = .data[[byvar]], yend = .data[[byvar]],
                     x = 1, xend = ratio, color = color))+
    ggtitle(title) +
    xlab("Ratio") +
    scale_color_manual(values = c("deepskyblue", "gray40")) +
    theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          panel.border = element_blank(),
          panel.spacing.y = unit(0.2, "line"),
          strip.background = element_rect(colour = NA, fill = "white"),
          strip.text = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = axis.text.y.size),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  return(list(g_out = g_out, plot_data = ratio_ypll_death))
}
