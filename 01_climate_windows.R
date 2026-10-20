# Part 1: Climate windows and Figure 1a-e
# Required packages: dplyr, ggplot2, climwin, lme4, terra, scales.
# Download WorldClim rasters separately and place them below worldclim_dir.

public_dir <- "."
data_dir <- file.path(public_dir, "data")
output_dir <- file.path(public_dir, "results")
figure_dir <- file.path(output_dir, "figures")
worldclim_dir <- file.path(public_dir, "external_data", "worldclim", "historical_tmin")

# -----------------------------------------------------------------------------
# Part 0: Create output folders.
# -----------------------------------------------------------------------------
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Part 1: Summarize laying dates and test differences among years.
# -----------------------------------------------------------------------------

library(dplyr)

dx_data <- read.csv(
  file.path(data_dir, "dxfamily.csv"),
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Keep records with a year and laying date.
lay_data <- dx_data %>%
  mutate(
    year = as.numeric(year),
    ld41 = as.numeric(ld41)
  ) %>%
  filter(
    !is.na(year),
    !is.na(ld41)
  )

cat("Breeding records: ", nrow(lay_data), "\n")
cat("Females: ", n_distinct(lay_data$female), "\n")
cat("Study years: ", min(lay_data$year), "-", max(lay_data$year), "\n")
cat("Number of years: ", n_distinct(lay_data$year), "\n\n")

# Calculate annual means, standard deviations, and sample sizes.
annual_laying <- lay_data %>%
  group_by(year) %>%
  summarise(
    n_records = n(),
    n_females = n_distinct(female),
    mean_ld41 = mean(ld41, na.rm = TRUE),
    sd_ld41 = sd(ld41, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year)

print(annual_laying)

# Convert days after 1 April to calendar dates for reporting.
origin_date <- as.Date("2001-04-01")

# Round to the nearest day.
round_day <- function(x) {
  floor(x + 0.5)
}

to_calendar_date <- function(x) {
  converted_date <- origin_date + round_day(x)
  paste(
    as.integer(format(converted_date, "%d")),
    format(converted_date, "%B")
  )
}

annual_laying$mean_calendar_date <- to_calendar_date(annual_laying$mean_ld41)

# Summarize the 15 annual means, with equal weight per year.
annual_mean <- mean(annual_laying$mean_ld41)
annual_sd <- sd(annual_laying$mean_ld41)
annual_min <- min(annual_laying$mean_ld41)
annual_max <- max(annual_laying$mean_ld41)

cat("Mean of annual means: ", round(annual_mean, 2), "days; ",
    to_calendar_date(annual_mean), "\n", sep = "")
cat("SD among annual means: ", round(annual_sd, 2), " days\n", sep = "")
cat("Range of annual means: ", round(annual_min, 2), "-", round(annual_max, 2),
    " days; ", to_calendar_date(annual_min), "-",
    to_calendar_date(annual_max), "\n", sep = "")

# One-way ANOVA of laying dates among years.
lay_data$year_factor <- factor(lay_data$year)

anova_model <- aov(ld41 ~ year_factor, data = lay_data)
anova_result <- summary(anova_model)[[1]]

df_year <- anova_result["year_factor", "Df"]
df_residual <- anova_result["Residuals", "Df"]
F_value <- anova_result["year_factor", "F value"]
P_value <- anova_result["year_factor", "Pr(>F)"]

print(anova_result)
cat("ANOVA: F(", df_year, ", ", df_residual, ") = ",
    round(F_value, 2), ", P = ", format.pval(P_value, digits = 4),
    "\n", sep = "")

# Calculate the 5th and 95th percentiles across breeding records.
record_quantile <- quantile(
  lay_data$ld41,
  probs = c(0.05, 0.95),
  na.rm = TRUE,
  type = 7
)

q05 <- as.numeric(record_quantile[1])
q95 <- as.numeric(record_quantile[2])

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

anova_table <- data.frame(
  Test = "One-way ANOVA among 15 years",
  df_year = df_year,
  df_residual = df_residual,
  F_value = F_value,
  P_value = P_value
)
print(anova_table)
descriptive_table <- data.frame(
  Statistic = c(
    "Mean of 15 annual means",
    "SD among 15 annual means",
    "Minimum annual mean",
    "Maximum annual mean",
    "5th percentile of individual records",
    "95th percentile of individual records"
  ),
  Value_days = c(
    annual_mean,
    annual_sd,
    annual_min,
    annual_max,
    q05,
    q95
  ),
  Calendar_date = c(
    to_calendar_date(annual_mean),
    NA,
    to_calendar_date(annual_min),
    to_calendar_date(annual_max),
    to_calendar_date(q05),
    to_calendar_date(q95)
  )
)
print(descriptive_table)

print(annual_laying)
write.csv(
  annual_laying,
  file.path(output_dir, "part1_annual_laying_date_by_year.csv"),
  row.names = FALSE
)

print(anova_table)
write.csv(
  anova_table,
  file.path(output_dir, "part1_laying_date_anova.csv"),
  row.names = FALSE
)

print(descriptive_table)
write.csv(
  descriptive_table,
  file.path(output_dir, "part1_laying_date_descriptive_statistics.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 2: Figure 1a, annual laying-date distributions.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source(file.path(public_dir, "plot_style_publication.R"))

fig1a_data <- read.csv(
  file.path(data_dir, "dxfamily.csv"),
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(fig1a_data) <- sub("^\xef\xbb\xbf", "", names(fig1a_data), useBytes = TRUE)
names(fig1a_data) <- trimws(names(fig1a_data))

fig1a_data <- fig1a_data %>%
  mutate(
    year = as.numeric(year),
    ld41 = as.numeric(ld41)
  ) %>%
  filter(
    !is.na(year),
    !is.na(ld41)
  )

fig1a_years <- sort(unique(fig1a_data$year))
fig1a_breaks <- c(2005, 2010, 2015, 2020)

fig1a_box_color <- "#4E79A7"

p_fig1a <- ggplot(fig1a_data, aes(x = year, y = ld41, group = year)) +
  geom_boxplot(
    color = fig1a_box_color,
    fill = pub_alpha(fig1a_box_color, 0.14),
    width = 0.56,
    linewidth = 0.32,
    outlier.shape = NA
  ) +
  scale_x_continuous(
    breaks = fig1a_breaks,
    limits = c(min(fig1a_years) - 1.00, 2020 + 0.80),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(30, 80, by = 10),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  coord_cartesian(ylim = c(24, 82), clip = "off") +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8)

print(p_fig1a)

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

save_pub_pair(
  fix_pub_panel_size(p_fig1a),
  file.path(figure_dir, "fig1a_laying_date_by_year"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigure 1a saved to results/figures.\n")

# -----------------------------------------------------------------------------
# Part 3: Compare nine candidate climate variables using climwin.
# -----------------------------------------------------------------------------

library(climwin)
library(lme4)
library(dplyr)

clim <- read.csv(
  file.path(data_dir, "dangxiong_monthly.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

biol <- read.csv(
  file.path(data_dir, "dxfamily.csv"),
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove surrounding whitespace from climate column names.
names(clim) <- trimws(names(clim))

biol$female <- as.factor(biol$female)
biol$year <- as.factor(biol$year)
biol$ld41 <- as.numeric(biol$ld41)

biol <- biol %>%
  filter(
    !is.na(female),
    !is.na(year),
    !is.na(ld41),
    !is.na(date)
  )

cat("========== climwin input checks ==========\n")
cat("Climate months: ", nrow(clim), "\n")
cat("Breeding records: ", nrow(biol), "\n")
cat("Number of years: ", n_distinct(biol$year), "\n")
cat("female class: ", class(biol$female), "\n")
cat("year class: ", class(biol$year), "\n")

# Candidate climate variables.
climate_variables <- data.frame(
  variable = c(
    "meantem",
    "mintem",
    "maxtem",
    "sumpre",
    "meanpre",
    "maxdaypre",
    "extmintem",
    "extmaxtem",
    "meanwin"
  ),
  label = c(
    "Mean temperature",
    "Mean minimum temperature",
    "Mean maximum temperature",
    "Accumulated precipitation",
    "Mean precipitation",
    "Maximum daily precipitation",
    "Extreme minimum temperature",
    "Extreme maximum temperature",
    "Mean wind speed"
  ),
  stringsAsFactors = FALSE
)

missing_climate_vars <- setdiff(climate_variables$variable, names(clim))

if (length(missing_climate_vars) > 0) {
  stop(paste("Missing climate variables: ", paste(missing_climate_vars, collapse = ", ")))
}

# Fit a baseline model with female and year random intercepts.
baseline_model <- lme4::lmer(
  ld41 ~ (1 | female) + (1 | year),
  data = biol
)

# Set the number of randomizations per climate variable.
climwin_random_repeats <- 100

# Fit and randomize one climate-window analysis.
run_one_climwin <- function(varname) {

  cat("\nAnalyzing climate variable: ", varname, "\n")

  sliding_result <- slidingwin(
    xvar = list(Temp = clim[[varname]]),
    exclude = NA,
    cdate = clim$date,
    bdate = biol$date,
    baseline = baseline_model,
    range = c(12, 0),
    type = "absolute",
    refday = c(1, 6),
    stat = "mean",
    func = "lin",
    cmissing = "method2",
    cinterval = "month"
  )

  random_result <- randwin(
    repeats = climwin_random_repeats,
    xvar = list(Temp = clim[[varname]]),
    exclude = NA,
    cdate = clim$date,
    bdate = biol$date,
    baseline = baseline_model,
    range = c(12, 0),
    type = "absolute",
    refday = c(1, 6),
    stat = "mean",
    func = "lin",
    cmissing = "method2",
    cinterval = "month"
  )

  list(
    sens = sliding_result,
    rand = random_result
  )
}

climwin_results_list <- lapply(
  climate_variables$variable,
  run_one_climwin
)

names(climwin_results_list) <- climate_variables$variable

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

saveRDS(
  climwin_results_list,
  file.path(output_dir, "part1_climwin_results.rds")
)

cat("\nclimwin results saved to results/part1_climwin_results.rds\n")

# -----------------------------------------------------------------------------
# Part 4: Summarize the best climate window for each variable.
# -----------------------------------------------------------------------------

library(climwin)
library(dplyr)

climwin_results_list <- readRDS(file.path(output_dir, "part1_climwin_results.rds"))

biol <- read.csv(
  file.path(data_dir, "dxfamily.csv"),
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

biol$year <- as.factor(biol$year)

climate_variables <- data.frame(
  variable = c(
    "meantem",
    "mintem",
    "maxtem",
    "sumpre",
    "meanpre",
    "maxdaypre",
    "extmintem",
    "extmaxtem",
    "meanwin"
  ),
  label = c(
    "Mean temperature",
    "Mean minimum temperature",
    "Mean maximum temperature",
    "Accumulated precipitation",
    "Mean precipitation",
    "Maximum daily precipitation",
    "Extreme minimum temperature",
    "Extreme maximum temperature",
    "Mean wind speed"
  ),
  stringsAsFactors = FALSE
)

# Convert window offsets into calendar-month labels.
window_to_month <- function(open, close, ref_month = 6) {

  mb_open <- abs(as.integer(open))
  mb_close <- abs(as.integer(close))

  m1 <- ((ref_month - mb_open - 1) %% 12) + 1
  m2 <- ((ref_month - mb_close - 1) %% 12) + 1

  paste0(m1, "-", m2, " (months)")
}

# Extract and format the randomization P value.
get_p_value <- function(dataset, datasetrand, sample_size) {

  p_raw <- tryCatch(
    pvalue(
      dataset = dataset,
      datasetrand = datasetrand,
      metric = "AIC",
      sample.size = sample_size
    ),
    error = function(e) NA
  )

  p_raw_text <- as.character(p_raw[1])

  if (is.na(p_raw_text)) {
    p_number <- NA_real_
    p_text <- NA_character_

  } else if (grepl("^<", p_raw_text)) {
    p_number <- as.numeric(gsub("<", "", p_raw_text))
    p_text <- p_raw_text

  } else {
    p_number <- suppressWarnings(as.numeric(p_raw_text))

    if (is.na(p_number)) {
      p_text <- NA_character_
    } else if (p_number < 0.001) {
      p_text <- "<0.001"
    } else if (p_number < 0.01) {
      p_text <- "<0.01"
    } else if (p_number < 0.05) {
      p_text <- "<0.05"
    } else {
      p_text <- sprintf("%.3f", p_number)
    }
  }

  list(
    p_text = p_text,
    p_number = p_number
  )
}

# Extract the minimum-deltaAICc window and its model statistics.
extract_one_climwin <- function(varname) {

  sens <- climwin_results_list[[varname]]$sens
  rand <- climwin_results_list[[varname]]$rand

  dataset <- sens[[1]]$Dataset
  best <- dataset[which.min(dataset$deltaAICc), , drop = FALSE]

  p_value <- get_p_value(
    dataset = dataset,
    datasetrand = rand[[1]],
    sample_size = n_distinct(biol$year)
  )

  data.frame(
    variable = varname,
    window_open = best$WindowOpen,
    window_close = best$WindowClose,
    window_month = window_to_month(
      best$WindowOpen,
      best$WindowClose,
      ref_month = best$Reference.month
    ),
    beta = best$ModelBeta,
    se = best$Std.Error,
    lower_95CI = best$ModelBeta - 1.96 * best$Std.Error,
    upper_95CI = best$ModelBeta + 1.96 * best$Std.Error,
    deltaAICc = best$deltaAICc,
    weight = best$ModWeight,
    p_PDAICc = p_value$p_text,
    p_PDAICc_number = p_value$p_number
  )
}

climwin_summary <- bind_rows(
  lapply(climate_variables$variable, extract_one_climwin)
) %>%
  left_join(climate_variables, by = "variable") %>%
  arrange(p_PDAICc_number, deltaAICc) %>%
  mutate(
    beta_se = sprintf("%.2f ± %.2f", beta, se),
    CI = sprintf("%.2f to %.2f", lower_95CI, upper_95CI)
  ) %>%
  select(
    variable,
    label,
    window_open,
    window_close,
    window_month,
    beta,
    se,
    beta_se,
    lower_95CI,
    upper_95CI,
    CI,
    deltaAICc,
    weight,
    p_PDAICc,
    p_PDAICc_number
  )

print(climwin_summary)

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

write.csv(
  climwin_summary,
  file.path(output_dir, "part1_climwin_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 5: Figure 1b, minimum-temperature climate-window heatmap.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source(file.path(public_dir, "plot_style_publication.R"))

# Find the Dataset element in the climwin sensitivity results.
extract_climwin_dataset <- function(climwin_results_list, climate_name = "mintem") {
  climate_result <- climwin_results_list[[climate_name]]

  if (is.null(climate_result) || is.null(climate_result[["sens"]])) {
    stop("Cannot find climwin sensitivity result for ", climate_name)
  }

  sens_items <- climate_result[["sens"]]

  dataset_index <- which(vapply(
    sens_items,
    function(x) is.list(x) && !is.null(x[["Dataset"]]),
    logical(1)
  ))

  if (!length(dataset_index)) {
    stop("Cannot find Dataset inside ", climate_name, " sensitivity result")
  }

  sens_items[[dataset_index[1]]][["Dataset"]]
}

climwin_results_list <- readRDS(file.path(output_dir, "part1_climwin_results.rds"))

main_climate_variable <- "mintem"

main_dataset <- extract_climwin_dataset(
  climwin_results_list,
  climate_name = main_climate_variable
)

# Identify the window with the lowest deltaAICc.
fig1b_data <- main_dataset %>%
  mutate(
    WindowOpen = as.numeric(WindowOpen),
    WindowClose = as.numeric(WindowClose),
    deltaAICc = as.numeric(deltaAICc)
  ) %>%
  filter(
    !is.na(WindowOpen),
    !is.na(WindowClose),
    !is.na(deltaAICc)
  )

best <- fig1b_data[which.min(fig1b_data$deltaAICc), , drop = FALSE]

p_fig1b <- ggplot(
  fig1b_data,
  aes(x = WindowClose, y = WindowOpen, fill = deltaAICc)
) +
  geom_tile(color = "white", linewidth = 0.06) +
  geom_point(
    data = best,
    aes(x = WindowClose, y = WindowOpen),
    inherit.aes = FALSE,
    shape = 23,
    size = 1.65,
    stroke = 0.34,
    fill = pub_palette[["purple"]],
    color = "white"
  ) +
  scale_fill_gradient2(
    low = pub_palette[["climate_red"]],
    mid = pub_palette[["very_light_grey"]],
    high = pub_palette[["climate_blue"]],
    midpoint = 0,
    breaks = c(-4, -2, 0, 2),
    labels = function(x) sprintf("%2.0f", x),
    name = expression(Delta*"AICc")
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(7.0, "mm"),
      barwidth = grid::unit(1.7, "mm"),
      ticks.linewidth = 0.20,
      label.hjust = 1
    )
  ) +
  theme(
    legend.position = c(0.92, 0.40),
    legend.justification = c(0.5, 0.5),
    legend.title = element_text(size = 4.9),
    legend.text = element_text(size = 5.0, hjust = 1),
    plot.margin = margin(2.5, 2.0, 2.5, 2.5)
  )

print(p_fig1b)

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

save_pub_pair(
  fix_pub_panel_size(p_fig1b),
  file.path(figure_dir, "fig1b_climate_window_mintem_deltaAICc"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigure 1b saved to results/figures.\n")

# -----------------------------------------------------------------------------
# Part 6: Figure 1c from the climwin best-window data.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source(file.path(public_dir, "plot_style_publication.R"))

fig1c_line_color <- "#4E79A7"
fig1c_ribbon_color <- "#D7E3EE"

# Find the BestModelData element in the climwin sensitivity results.
extract_climwin_best_model_data <- function(climwin_results_list, climate_name = "mintem") {
  climate_result <- climwin_results_list[[climate_name]]

  if (is.null(climate_result) || is.null(climate_result[["sens"]])) {
    stop("Cannot find climwin sensitivity result for ", climate_name)
  }

  sens_items <- climate_result[["sens"]]

  best_model_index <- which(vapply(
    sens_items,
    function(x) is.list(x) && !is.null(x[["BestModelData"]]),
    logical(1)
  ))

  if (!length(best_model_index)) {
    stop("Cannot find BestModelData inside ", climate_name, " sensitivity result")
  }

  sens_items[[best_model_index[1]]][["BestModelData"]]
}

climwin_results_list <- readRDS(file.path(output_dir, "part1_climwin_results.rds"))

main_climate_variable <- "mintem"

main_best_model_data <- extract_climwin_best_model_data(
  climwin_results_list,
  climate_name = main_climate_variable
)

fig1c_data <- main_best_model_data %>%
  mutate(
    year = as.numeric(as.character(year)),
    yvar = as.numeric(yvar),
    climate = as.numeric(climate)
  ) %>%
  filter(
    !is.na(year),
    !is.na(yvar),
    !is.na(climate)
  ) %>%
  group_by(year) %>%
  summarise(
    mean_laying_date = mean(yvar, na.rm = TRUE),
    winter_temperature = mean(climate, na.rm = TRUE),
    n_records = n(),
    .groups = "drop"
  ) %>%
  arrange(year)

cat("\n========== Figure 1c plotting data ==========\n")
print(fig1c_data)

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

write.csv(
  fig1c_data,
  file.path(output_dir, "part1_fig1c_laying_date_winter_temperature_data.csv"),
  row.names = FALSE
)

p_fig1c <- ggplot(
  fig1c_data,
  aes(x = winter_temperature, y = mean_laying_date)
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = pub_scatter$line_width,
    color = fig1c_line_color,
    fill = pub_alpha(fig1c_ribbon_color, 0.75),
    alpha = 0.6
  ) +
  geom_point(
    size = pub_scatter$point_size,
    stroke = pub_scatter$point_stroke,
    shape = 21,
    fill = "white",
    color = fig1c_line_color,
    alpha = 0.88
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8)

print(p_fig1c)

save_pub_pair(
  fix_pub_panel_size(p_fig1c),
  file.path(figure_dir, "fig1c_laying_date_winter_mintem"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigure 1c saved from climwin best-window data.\n")

# -----------------------------------------------------------------------------
# Part 7: Figure 1e, annual mean laying-date trend, 2005-2019.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source(file.path(public_dir, "plot_style_publication.R"))

fig1e_line_color <- "#4E79A7"
fig1e_ribbon_color <- "#D7E3EE"

fig1e_data <- read.csv(
  file.path(data_dir, "dxfamily.csv"),
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fig1e_data <- fig1e_data %>%
  mutate(
    year = as.numeric(year),
    ld41 = as.numeric(ld41)
  ) %>%
  filter(
    !is.na(year),
    !is.na(ld41),
    year >= 2005,
    year <= 2019
  )

fig1e_annual <- fig1e_data %>%
  group_by(year) %>%
  summarise(
    mean_laying_date = mean(ld41, na.rm = TRUE),
    se_laying_date = sd(ld41, na.rm = TRUE) / sqrt(n()),
    n_records = n(),
    .groups = "drop"
  ) %>%
  arrange(year)

cat("\n========== Figure 1e plotting data ==========\n")
print(fig1e_annual)

# Fit a linear trend to the 15 annual means.
fig1e_lm <- lm(
  mean_laying_date ~ year,
  data = fig1e_annual
)

print(fig1e_lm)

fig1e_lm_summary <- summary(fig1e_lm)

fig1e_slope <- coef(fig1e_lm)[["year"]]
fig1e_p <- fig1e_lm_summary$coefficients["year", "Pr(>|t|)"]
fig1e_r <- cor(fig1e_annual$year, fig1e_annual$mean_laying_date)

cat("\n========== Figure 1e annual mean laying-date trend ==========\n")
cat("r = ", round(fig1e_r, 3), "\n", sep = "")
cat("slope = ", round(fig1e_slope, 3), " days/year\n", sep = "")
cat("P = ", format.pval(fig1e_p, digits = 4), "\n", sep = "")

figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

write.csv(
  fig1e_annual,
  file.path(output_dir, "part1_fig1e_annual_laying_date_trend_data.csv"),
  row.names = FALSE
)

fig1e_trend_table <- data.frame(
  Test = "Temporal trend in annual mean laying date",
  r = fig1e_r,
  slope_days_per_year = fig1e_slope,
  P_value = fig1e_p
)

print(fig1e_trend_table)

# Export the trend estimates.
write.csv(
  fig1e_trend_table,
  file.path(output_dir, "part1_fig1e_annual_laying_date_trend_stats.csv"),
  row.names = FALSE
)

p_fig1e <- ggplot(
  fig1e_annual,
  aes(x = year, y = mean_laying_date)
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = pub_scatter$line_width,
    linetype = "dashed",
    color = fig1e_line_color,
    fill = pub_alpha(fig1e_ribbon_color, 0.75),
    alpha = 0.6
  ) +
  geom_point(
    shape = 21,
    size = pub_scatter$point_size,
    stroke = pub_scatter$point_stroke,
    color = fig1e_line_color,
    fill = "white",
    alpha = 0.88
  ) +
  scale_x_continuous(
    breaks = c(2005, 2010, 2015, 2020),
    limits = c(2004, 2020.8),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8) +
  theme(
    legend.position = "none"
  )

print(p_fig1e)

save_pub_pair(
  fix_pub_panel_size(p_fig1e),
  file.path(figure_dir, "fig1e_annual_laying_date_trend"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigure 1e saved to results/figures.\n")

# -----------------------------------------------------------------------------
# Part 8: Figure 1d and temporal trend summaries.
# -----------------------------------------------------------------------------

# Winter for breeding year Y is the mean of November and December in Y-1 and January in Y.
library(terra)
library(dplyr)
library(ggplot2)

source(file.path(public_dir, "plot_style_publication.R"))

fig1d_line_color <- "#4E79A7"
fig1d_ribbon_color <- "#D7E3EE"
fig1d_marker_color <- "#D99045"

# Define the study site (WGS84).
site <- terra::vect(
  data.frame(lon = 91.09861667, lat = 30.49945),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

# Extract WorldClim monthly temperatures and save them for reuse (first run).
worldclim_files <- list.files(
  worldclim_dir,
  pattern = "\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(worldclim_files) == 0) {
  stop("No WorldClim GeoTIFFs found. Check worldclim_dir in the configuration section.")
}

# Read year and month from the raster filenames.
worldclim_file_info <- data.frame(
  file = worldclim_files,
  fname = basename(worldclim_files),
  stringsAsFactors = FALSE
) %>%
  mutate(
    year = as.numeric(sub(".*tmin_([0-9]{4})-([0-9]{2})\\.tif$", "\\1", fname)),
    month = as.numeric(sub(".*tmin_([0-9]{4})-([0-9]{2})\\.tif$", "\\2", fname))
  ) %>%
  filter(
    !is.na(year),
    !is.na(month)
  ) %>%
  arrange(year, month)

cat("GeoTIFF files: ", nrow(worldclim_file_info), "\n")
cat("Year range: ", min(worldclim_file_info$year), "-", max(worldclim_file_info$year), "\n")

# Select the months required for winters in 1960-2019.
worldclim_needed <- worldclim_file_info %>%
  filter(
    (year >= 1959 & year <= 2018 & month %in% c(11, 12)) |
      (year >= 1960 & year <= 2019 & month == 1)
  )

cat("Monthly rasters selected: ", nrow(worldclim_needed), "\n")

# Extract monthly minimum temperatures at the study site.
worldclim_needed$tmin <- sapply(
  worldclim_needed$file,
  function(f) {
    r <- terra::rast(f)
    terra::extract(r, site)[1, 2]
  }
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save extracted monthly values and source metadata.
worldclim_site_tmin_needed <- worldclim_needed %>%
  select(year, month, fname, tmin) %>%
  mutate(
    lon = 91.09861667,
    lat = 30.49945,
    source = "WorldClim historical monthly tmin",
    note = "Months needed for Nov-Dec-Jan winter minimum temperature"
  )

saveRDS(
  worldclim_site_tmin_needed,
  file.path(output_dir, "part1_worldclim_site_tmin_needed_months_1959_2019.rds")
)

write.csv(
  worldclim_site_tmin_needed,
  file.path(output_dir, "part1_worldclim_site_tmin_needed_months_1959_2019.csv"),
  row.names = FALSE
)

cat("min = ", min(worldclim_needed$tmin, na.rm = TRUE), "\n", sep = "")
cat("max = ", max(worldclim_needed$tmin, na.rm = TRUE), "\n", sep = "")

# Read saved monthly temperatures; subsequent runs can resume here.
worldclim_needed <- readRDS(
  file.path(output_dir, "part1_worldclim_site_tmin_needed_months_1959_2019.rds")
)

# Retrieve one monthly value for the winter calculation.
get_tmin <- function(y, m) {

  value <- worldclim_needed$tmin[
    worldclim_needed$year == y &
      worldclim_needed$month == m
  ]

  if (length(value) == 0) {
    return(NA_real_)
  } else {
    return(value[1])
  }
}

# Average the three months without removing missing values.

winter_tmin_annual <- data.frame(
  year = 1960:2019
) %>%
  mutate(
    nov_previous_year = sapply(year, function(y) get_tmin(y - 1, 11)),
    dec_previous_year = sapply(year, function(y) get_tmin(y - 1, 12)),
    jan_current_year = sapply(year, function(y) get_tmin(y, 1)),
    winter_tmin = rowMeans(
      data.frame(nov_previous_year, dec_previous_year, jan_current_year),
      na.rm = FALSE
    ),
    study_period = ifelse(year >= 2005 & year <= 2019, "2005-2019", "1960-2004")
  )

print(head(winter_tmin_annual))
print(tail(winter_tmin_annual))
cat("Missing winter_tmin values: ", sum(is.na(winter_tmin_annual$winter_tmin)), "\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

print(winter_tmin_annual)
write.csv(
  winter_tmin_annual,
  file.path(output_dir, "part1_fig1d_worldclim_winter_tmin_1960_2019.csv"),
  row.names = FALSE
)

# Test the long-term winter-temperature trend using WorldClim.
long_temp_test <- cor.test(
  winter_tmin_annual$year,
  winter_tmin_annual$winter_tmin,
  method = "pearson"
)

long_temp_r <- as.numeric(long_temp_test$estimate)
long_temp_p <- long_temp_test$p.value

# Use the local monthly climate data for the 2005-2019 winter-temperature trend.
local_clim_for_study <- read.csv(
  file.path(data_dir, "dangxiong_monthly.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(local_clim_for_study) <- trimws(names(local_clim_for_study))

get_local_mintem <- function(y, m) {

  value <- local_clim_for_study$mintem[
    local_clim_for_study$year == y &
      local_clim_for_study$month == m
  ]

  if (length(value) == 0) {
    return(NA_real_)
  } else {
    return(value[1])
  }
}

study_temp_data <- data.frame(
  year = 2005:2019
) %>%
  mutate(
    nov_previous_year = sapply(year, function(y) get_local_mintem(y - 1, 11)),
    dec_previous_year = sapply(year, function(y) get_local_mintem(y - 1, 12)),
    jan_current_year = sapply(year, function(y) get_local_mintem(y, 1)),
    winter_tmin = rowMeans(
      data.frame(nov_previous_year, dec_previous_year, jan_current_year),
      na.rm = FALSE
    )
  )

head(study_temp_data)

write.csv(
  study_temp_data,
  file.path(output_dir, "part1_study_period_local_winter_tmin_2005_2019.csv"),
  row.names = FALSE
)

study_temp_test <- cor.test(
  study_temp_data$year,
  study_temp_data$winter_tmin,
  method = "pearson"
)

study_temp_r <- as.numeric(study_temp_test$estimate)
study_temp_p <- study_temp_test$p.value

# Summarize temperature trends and reuse the Figure 1e trend results.
trend_results <- data.frame(
  analysis = c(
    "Winter minimum temperature trend, 1960-2019",
    "Winter minimum temperature trend, 2005-2019",
    "Annual mean breeding time trend, 2005-2019"
  ),
  r = c(
    long_temp_r,
    study_temp_r,
    fig1e_r
  ),
  P_value = c(
    long_temp_p,
    study_temp_p,
    fig1e_p
  )
)

print(trend_results)

# Plot the WorldClim winter-temperature series (1960-2019).
figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

fig1d_y_range <- range(winter_tmin_annual$winter_tmin, na.rm = TRUE)
fig1d_y_span <- diff(fig1d_y_range)

fig1d_study_y_range <- range(
  winter_tmin_annual$winter_tmin[
    winter_tmin_annual$year >= 2005 &
      winter_tmin_annual$year <= 2019
  ],
  na.rm = TRUE
)

fig1d_marker <- data.frame(
  year = c(2005, 2019),
  y = fig1d_study_y_range[1] - 0.03 * fig1d_y_span,
  yend = fig1d_study_y_range[2] + 0.03 * fig1d_y_span
)

p_fig1d <- ggplot(
  winter_tmin_annual,
  aes(x = year, y = winter_tmin)
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = fig1d_line_color,
    linewidth = pub_scatter$line_width,
    fill = pub_alpha(fig1d_ribbon_color, 0.75),
    alpha = 0.6

  ) +
  geom_point(
    shape = 21,
    size = pub_scatter$point_size,
    stroke = pub_scatter$point_stroke,
    color = fig1d_line_color,
    fill = "white",
    alpha = 0.88
  ) +
  scale_x_continuous(
    breaks = c(1960, 1980, 2000, 2020),
    limits = c(1957, 2025),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8) +
  theme(
    legend.position = "none"
  )

print(p_fig1d)

save_pub_pair(
  fix_pub_panel_size(p_fig1d),
  file.path(figure_dir, "fig1d_long_term_winter_tmin_trend"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigure 1d saved to results/figures.\n")
