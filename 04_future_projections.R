# Part 4: Future warming, adaptive capacity, and Figure 3
# Required packages: terra, dplyr, stringr, purrr, coda, ggplot2, tidyr, cowplot, scales.
# Download WorldClim baseline and future minimum-temperature rasters separately.

library(terra)
library(dplyr)
library(stringr)
library(purrr)
library(coda)
library(ggplot2)
library(tidyr)

# -----------------------------------------------------------------------------
# Part 0: Set paths and fixed inputs.
# -----------------------------------------------------------------------------
public_dir <- "."
data_dir <- file.path(public_dir, "data")
output_dir <- file.path(public_dir, "results")
figure_dir <- file.path(output_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

current_dir <- file.path(public_dir, "external_data", "worldclim", "current_tmin")
future_dir <- file.path(public_dir, "external_data", "worldclim", "future_tmin")

data_file <- file.path(data_dir, "dxfamily.csv")
model_file <- file.path(output_dir, "mod_GE.rds")
part3_R_rds_file <- file.path(output_dir, "part3_evolutionary_response_posterior_samples.rds")
part3_R_csv_file <- file.path(output_dir, "part3_evolutionary_response_posterior_samples.csv")

# Start accumulating evolutionary change from 2020.
baseline_year <- 2020

# Female generation time in years per generation.
gen_time <- 1.89

# Use up to 3,000 uncertainty draws.
n_draw_wanted <- 3000

scenario_levels <- c("ssp126", "ssp245", "ssp370", "ssp585")
scenario_labels <- c(
  "ssp126" = "SSP1-2.6",
  "ssp245" = "SSP2-4.5",
  "ssp370" = "SSP3-7.0",
  "ssp585" = "SSP5-8.5"
)

# Summarize draws with quantiles and shortest 95% intervals.
q025 <- function(x) as.numeric(quantile(x, 0.025, na.rm = TRUE))
q975 <- function(x) as.numeric(quantile(x, 0.975, na.rm = TRUE))

posterior_summary <- function(x, parameter) {
  h <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    Parameter = parameter,
    Mean = mean(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    HPD_lower = as.numeric(h[1]),
    HPD_upper = as.numeric(h[2]),
    q025 = q025(x),
    q975 = q975(x),
    Eff.samp = as.numeric(effectiveSize(as.mcmc(x))),
    row.names = NULL
  )
}

# Define the study site (WGS84).
site <- vect(
  data.frame(lon = 91.09861667, lat = 30.49945),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

# -----------------------------------------------------------------------------
# Part 1: Extract the WorldClim baseline winter temperature.
# -----------------------------------------------------------------------------

# Read the 12 baseline monthly minimum-temperature rasters.
current_files <- list.files(
  current_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)
current_meta <- data.frame(file = current_files) %>%
  mutate(
    fname = basename(file),
    month = str_extract(fname, "(?<=tmin_)[0-9]{2}") %>% as.numeric()
  ) %>%
  filter(!is.na(month)) %>%
  arrange(month)

if (nrow(current_meta) < 12) {
  stop(paste0(
    "Fewer than 12 baseline climate files found. Check: ",
    current_dir
  ))
}

# Extract monthly minimum temperatures at the study site.
current_meta$value <- map_dbl(current_meta$file, function(f) {
  r <- rast(f)
  val <- terra::extract(r, site)[1, 2]
  as.numeric(val)
})

# Average January, November, and December for the 1970-2000 baseline.
current_tem11_1 <- current_meta %>%
  filter(month %in% c(1, 11, 12)) %>%
  summarise(tem11_1 = mean(value, na.rm = TRUE)) %>%
  pull(tem11_1)

cat("\nBaseline winter temperature = ", round(current_tem11_1, 3), " degrees Celsius\n")

current_temperature_output <- current_meta %>%
  mutate(
    current_tem11_1 = current_tem11_1,
    source = "WorldClim current baseline 1970-2000"
  )
print(current_temperature_output)
write.csv(
  current_temperature_output,
  file.path(output_dir, "part4_worldclim_current_temperature.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 2: Extract future winter temperatures and warming increments.
# -----------------------------------------------------------------------------

# Read future files containing 12 monthly layers per GCM, SSP scenario, and period.
future_files <- list.files(
  future_dir,
  pattern = "^wc2\\.1_30s_tmin_.*_ssp(126|245|370|585)_[0-9]{4}-[0-9]{4}\\.tif$",
  full.names = TRUE
)

if (length(future_files) == 0) {
  stop(paste0(
    "No future climate files found. Check: ",
    future_dir
  ))
}

# Extract GCM, scenario, and period identifiers from filenames.
future_meta <- data.frame(file = future_files) %>%
  mutate(
    fname = basename(file),
    scenario = str_extract(fname, "ssp126|ssp245|ssp370|ssp585"),
    scenario_lab = recode(scenario, !!!scenario_labels),
    period = str_extract(fname, "[0-9]{4}-[0-9]{4}"),
    period_start = as.numeric(str_sub(period, 1, 4)),
    period_end = as.numeric(str_sub(period, 6, 9)),
    period_mid = (period_start + period_end) / 2,
    GCM = fname %>%
      str_remove("^wc2\\.1_30s_tmin_") %>%
      str_remove("_ssp(126|245|370|585)_[0-9]{4}-[0-9]{4}\\.tif$")
  ) %>%
  arrange(scenario, period_mid, GCM)

# Calculate winter temperature from monthly layers 1, 11, and 12.
extract_tem11_1_from_tn <- function(f, site) {
  r <- rast(f)

  if (nlyr(r) < 12) {
    stop(paste("Fewer than 12 monthly layers in: ", f))
  }

  vals <- terra::extract(r, site)[1, -1]
  vals <- as.numeric(vals)

  mean(vals[c(1, 11, 12)], na.rm = TRUE)
}

future_meta$future_tem11_1 <- map_dbl(
  future_meta$file,
  extract_tem11_1_from_tn,
  site = site
)

# Calculate warming relative to the WorldClim baseline.
clim_extract_df <- future_meta %>%
  mutate(
    current_tem11_1 = current_tem11_1,
    dT = future_tem11_1 - current_tem11_1,
    scenario = factor(scenario, levels = scenario_levels),
    scenario_lab = factor(
      scenario_lab,
      levels = unname(scenario_labels[scenario_levels])
    )
  ) %>%
  select(
    GCM,
    scenario,
    scenario_lab,
    period,
    period_start,
    period_end,
    period_mid,
    current_tem11_1,
    future_tem11_1,
    dT,
    file
  )

write.csv(
  clim_extract_df,
  file.path(output_dir, "part4_worldclim_future_temperature_extract.csv"),
  row.names = FALSE
)

# Summarize warming across GCMs using the mean, minimum, and maximum.
clim_summary <- clim_extract_df %>%
  group_by(scenario, scenario_lab, period, period_mid, period_end) %>%
  summarise(
    mean_dT = mean(dT, na.rm = TRUE),
    min_dT = min(dT, na.rm = TRUE),
    max_dT = max(dT, na.rm = TRUE),
    n_GCM = n(),
    .groups = "drop"
  ) %>%
  arrange(scenario, period_mid)

write.csv(
  clim_summary,
  file.path(output_dir, "part4_worldclim_future_temperature_summary.csv"),
  row.names = FALSE
)

# Save the extracted temperatures for reuse.
clim_cache_file <- file.path(output_dir, "part4_worldclim_temperature_cache.rds")

saveRDS(
  list(
    current_meta = current_meta,
    current_tem11_1 = current_tem11_1,
    current_temperature_output = current_temperature_output,
    future_meta = future_meta,
    clim_extract_df = clim_extract_df,
    clim_summary = clim_summary,
    site = site,
    current_dir = current_dir,
    future_dir = future_dir,
    created_at = Sys.time()
  ),
  clim_cache_file,
  compress = "xz"
)

# Read saved climate extractions; subsequent runs can resume here.
clim_cache_file <- file.path(output_dir, "part4_worldclim_temperature_cache.rds")
clim_cache <- readRDS(clim_cache_file)

current_meta <- clim_cache$current_meta
current_tem11_1 <- clim_cache$current_tem11_1
current_temperature_output <- clim_cache$current_temperature_output
future_meta <- clim_cache$future_meta
clim_extract_df <- clim_cache$clim_extract_df
clim_summary <- clim_cache$clim_summary

# -----------------------------------------------------------------------------
# Part 3: Calculate plastic, evolutionary, and combined capacities.
# -----------------------------------------------------------------------------

# Read breeding data and the random-regression animal model.
dat <- read.csv(
  data_file,
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mod_GE <- readRDS(model_file)

if (!"tem11_1" %in% names(dat)) {
  stop("Breeding data must contain tem11_1.")
}

if (!"tem11_1_z" %in% colnames(mod_GE$Sol)) {
  stop("The model must contain the fixed temperature effect tem11_1_z.")
}

dat$tem11_1 <- as.numeric(dat$tem11_1)

# Use the magnitude of the population temperature slope (days per degree Celsius).
beta_temp_post <- as.numeric(mod_GE$Sol[, "tem11_1_z"])
beta_abs_post <- abs(beta_temp_post)

# Calculate the historical temperature range across breeding records.
T_low <- q025(dat$tem11_1)
T_high <- q975(dat$tem11_1)
dT_hist <- T_high - T_low

cat("\nHistorical winter-temperature range (2.5th-97.5th percentiles) = ", round(dT_hist, 3), " degrees Celsius\n")

# Read the propagated heritability, selection, and response draws from Part 3.
if (file.exists(part3_R_rds_file)) {
  part3_R_post <- readRDS(part3_R_rds_file)
} else if (file.exists(part3_R_csv_file)) {
  part3_R_post <- read.csv(
    part3_R_csv_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  stop(paste0(
    "Part 3 evolutionary-response draw file not found: \n",
    part3_R_rds_file, "\n",
    part3_R_csv_file, "\n",
    "Run Part 3 to generate evolutionary-response draws."
  ))
}

if (!all(c("h2_int", "S", "R") %in% names(part3_R_post))) {
  stop("Part 3 draw files must contain h2_int, S, and R.")
}

part3_R_post <- part3_R_post %>%
  mutate(
    h2_int = as.numeric(h2_int),
    S = as.numeric(S),
    R = as.numeric(R)
  ) %>%
  filter(!is.na(h2_int), !is.na(S), !is.na(R))

if (nrow(part3_R_post) < 2) {
  stop("Too few usable Part 3 draws. Check the input file.")
}

# Sample temperature slopes and Part 3 draws to equal lengths.
set.seed(123)
n_draw <- min(n_draw_wanted, length(beta_abs_post), nrow(part3_R_post))

beta_abs_sample <- sample(beta_abs_post, n_draw, replace = FALSE)

part3_rows <- sample(seq_len(nrow(part3_R_post)), n_draw, replace = FALSE)
h2_sample <- part3_R_post$h2_int[part3_rows]
S_sample <- part3_R_post$S[part3_rows]
R_post <- part3_R_post$R[part3_rows]
R_abs_post <- abs(R_post)

# Calculate plastic capacity and evolutionary capacity accumulated by 2100.
plast_post <- beta_abs_sample * dT_hist

nyears_2100 <- 2100 - baseline_year
E_2100_post <- R_abs_post * nyears_2100 / gen_time
adapt_2100_post <- plast_post + E_2100_post

capacity_samples <- data.frame(
  beta_abs = beta_abs_sample,
  h2_int = h2_sample,
  S = S_sample,
  R = R_post,
  R_abs = R_abs_post,
  Dayplasticity = plast_post,
  Dayevolution_2100 = E_2100_post,
  Daytotal_2100 = adapt_2100_post
)
write.csv(
  capacity_samples,
  file.path(output_dir, "part4_capacity_posterior_samples.csv"),
  row.names = FALSE
)

# Summarize biological parameter and capacity draws.
capacity_summary <- bind_rows(
  posterior_summary(beta_temp_post, "Temperature slope beta_T (days per deg C)"),
  posterior_summary(beta_abs_post, "Absolute temperature slope abs(beta_T)"),
  posterior_summary(plast_post, "Plastic capacity"),
  posterior_summary(h2_sample, "Intercept heritability"),
  posterior_summary(S_sample, "Selection differential S"),
  posterior_summary(R_post, "Per-generation evolutionary response R"),
  posterior_summary(E_2100_post, "Evolutionary capacity by 2100"),
  posterior_summary(adapt_2100_post, "Total adaptive capacity by 2100")
)

capacity_summary$History_delta_T <- NA_real_
capacity_summary$History_delta_T[
  capacity_summary$Parameter == "Plastic capacity"
] <- dT_hist
cat("\n========== Adaptive-capacity summary ==========\n")
cat("Uncertainty draws: ", n_draw, "\n")
print(capacity_summary)
write.csv(
  capacity_summary,
  file.path(output_dir, "part4_adaptive_capacity_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 4: Summarize required advancement and matching probabilities.
# -----------------------------------------------------------------------------

# Calculate requirements, capacities, and matching probabilities separately for each GCM.
req_list <- lapply(seq_len(nrow(clim_extract_df)), function(i) {
  dT_i <- clim_extract_df$dT[i]
  period_end_i <- clim_extract_df$period_end[i]
  nyears_i <- period_end_i - baseline_year

  if (nyears_i < 0) {
    stop("A future period ends before baseline_year. Check the year inputs.")
  }

  req_post_i <- beta_abs_sample * dT_i

  E_post_i <- R_abs_post * nyears_i / gen_time

  adapt_post_i <- plast_post + E_post_i

  data.frame(
    GCM = clim_extract_df$GCM[i],
    scenario = as.character(clim_extract_df$scenario[i]),
    scenario_lab = as.character(clim_extract_df$scenario_lab[i]),
    period = clim_extract_df$period[i],
    period_mid = clim_extract_df$period_mid[i],
    period_end = period_end_i,
    dT = dT_i,

    mean_req = mean(req_post_i, na.rm = TRUE),
    lwr_req = q025(req_post_i),
    upr_req = q975(req_post_i),

    mean_plasticity = mean(plast_post, na.rm = TRUE),
    mean_evolution = mean(E_post_i, na.rm = TRUE),
    mean_total = mean(adapt_post_i, na.rm = TRUE),

    P_match = mean(adapt_post_i >= req_post_i, na.rm = TRUE),
    P_plasticity_match = mean(plast_post >= req_post_i, na.rm = TRUE),
    P_evolution_match = mean(E_post_i >= req_post_i, na.rm = TRUE)
  )
})
req_gcm_summary <- bind_rows(req_list) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    scenario_lab = factor(
      scenario_lab,
      levels = unname(scenario_labels[scenario_levels])
    )
  ) %>%
  arrange(scenario, period_mid, GCM)

write.csv(
  req_gcm_summary,
  file.path(output_dir, "part4_required_advance_by_gcm_period.csv"),
  row.names = FALSE
)

# Average GCM-specific matching probabilities within each scenario and period.
req_period_summary <- req_gcm_summary %>%
  rename(
    gcm_mean_req = mean_req,
    gcm_mean_plasticity = mean_plasticity,
    gcm_mean_evolution = mean_evolution,
    gcm_mean_total = mean_total
  ) %>%
  group_by(scenario, scenario_lab, period, period_mid, period_end) %>%
  summarise(
    n_GCM = n(),

    mean_dT = mean(dT, na.rm = TRUE),
    min_dT = min(dT, na.rm = TRUE),
    max_dT = max(dT, na.rm = TRUE),

    mean_req = mean(gcm_mean_req, na.rm = TRUE),
    min_req = min(gcm_mean_req, na.rm = TRUE),
    max_req = max(gcm_mean_req, na.rm = TRUE),

    mean_plasticity = mean(gcm_mean_plasticity, na.rm = TRUE),
    mean_evolution = mean(gcm_mean_evolution, na.rm = TRUE),
    mean_total = mean(gcm_mean_total, na.rm = TRUE),

    mean_P_match = mean(P_match, na.rm = TRUE),
    P_match_min = min(P_match, na.rm = TRUE),
    P_match_max = max(P_match, na.rm = TRUE),

    mean_P_plasticity_match = mean(P_plasticity_match, na.rm = TRUE),
    P_plasticity_match_min = min(P_plasticity_match, na.rm = TRUE),
    P_plasticity_match_max = max(P_plasticity_match, na.rm = TRUE),

    mean_P_evolution_match = mean(P_evolution_match, na.rm = TRUE),
    P_evolution_match_min = min(P_evolution_match, na.rm = TRUE),
    P_evolution_match_max = max(P_evolution_match, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(scenario, period_mid)

write.csv(
  req_period_summary,
  file.path(output_dir, "part4_required_advance_by_scenario_period.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 5: Figure 3, plasticity, evolution, and combined capacity.
# -----------------------------------------------------------------------------
source(file.path(public_dir, "plot_style_publication.R"))

library(cowplot)

# Set Figure 3 colors.
fig3_palette <- c(
  main_blue = "#4E79A7",
  ribbon_blue = "#D7E3EE",
  muted_purple = "#8E77B8",
  muted_teal = "#59A89C",
  muted_orange = "#D99045",
  phen_blue = "#4CA8BF",
  black = "#222222"
)

# Plot periods at their midpoints and label them by their end years.
period_axis <- req_period_summary %>%
  distinct(period_mid, period, period_end) %>%
  arrange(period_mid)

period_breaks <- period_axis$period_mid
period_labels_short <- as.character(period_axis$period_end)

fig3_x_limits <- c(min(period_breaks) - 7.0, max(period_breaks) + 7.0)

fig3_cols <- c(
  "Required" = fig3_palette[["black"]],
  "Plasticity" = fig3_palette[["main_blue"]],
  "Evolution" = fig3_palette[["muted_teal"]],
  "Total capacity" = fig3_palette[["muted_purple"]]
)

fig3_panel_levels <- c(
  "Plasticity",
  "Evolution",
  "Total capacity"
)

fig3_scenario_levels <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")

# Pool uncertainty draws across GCMs; report means and 2.5th-97.5th percentiles.
make_figure3_draw_summary <- function(rows) {
  draws <- bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    dT_i <- rows$dT[i]
    nyears_i <- rows$period_end[i] - baseline_year

    req_i <- beta_abs_sample * dT_i
    evo_i <- R_abs_post * nyears_i / gen_time
    plastic_i <- plast_post
    total_i <- plastic_i + evo_i

    data.frame(
      Dayrequired = req_i,
      Daytotal = total_i,
      Dayplasticity = plastic_i,
      Dayevolution = evo_i
    )
  }))

  draws %>%
    summarise(
      across(
        everything(),
        list(mean = ~ mean(.x, na.rm = TRUE), lwr = q025, upr = q975),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = c("component", ".value"),
      names_pattern = "(.+)_(mean|lwr|upr)"
    )
}

figure3_interval_data <- req_gcm_summary %>%
  group_by(scenario, scenario_lab, period, period_mid, period_end) %>%
  group_modify(~ make_figure3_draw_summary(.x), .keep = TRUE) %>%
  ungroup() %>%
  mutate(
    component_lab = recode(
      component,
      "Dayrequired" = "Required",
      "Daytotal" = "Total capacity",
      "Dayplasticity" = "Plasticity",
      "Dayevolution" = "Evolution"
    )
  )

# Compare each capacity with required advancement.
figure3_day_data <- bind_rows(
  figure3_interval_data %>%
    filter(component_lab %in% c("Required", "Plasticity")) %>%
    mutate(row_panel = "Plasticity"),
  figure3_interval_data %>%
    filter(component_lab %in% c("Required", "Evolution")) %>%
    mutate(row_panel = "Evolution"),
  figure3_interval_data %>%
    filter(component_lab %in% c("Required", "Total capacity")) %>%
    mutate(row_panel = "Total capacity")
) %>%
  transmute(
    scenario_lab,
    period_mid,
    period_end,
    row_panel,
    component_lab,
    mean,
    lwr,
    upr
  )

fig3_year_dodge_step <- 2.10

# Offset paired components horizontally within each future period.
figure3_plot_data <- figure3_day_data %>%
  mutate(
    row_panel = factor(row_panel, levels = fig3_panel_levels),
    scenario_lab = factor(scenario_lab, levels = fig3_scenario_levels),
    component_lab = factor(
      component_lab,
      levels = c("Required", "Plasticity", "Evolution", "Total capacity")
    )
  ) %>%
  group_by(row_panel, scenario_lab, period_mid) %>%
  mutate(
    component_n = n_distinct(component_lab),
    component_rank = match(
      as.numeric(component_lab),
      sort(unique(as.numeric(component_lab)))
    ),
    year_plot = period_mid + if_else(
      component_n > 1,
      (component_rank - (component_n + 1) / 2) * fig3_year_dodge_step,
      0
    )
  ) %>%
  ungroup()

write.csv(
  figure3_plot_data,
  file.path(output_dir, "part4_figure3_plot_data.csv"),
  row.names = FALSE
)

label_y_nudge <- 0.055

fig3_pmatch_label_size <- 2.40

figure3_panel_range <- figure3_plot_data %>%
  mutate(
    scenario_lab = as.character(scenario_lab),
    row_panel = as.character(row_panel)
  ) %>%
  group_by(scenario_lab, row_panel) %>%
  summarise(
    panel_min = min(c(lwr, mean, 0), na.rm = TRUE),
    panel_max = max(c(upr, mean, 0), na.rm = TRUE),
    panel_span = panel_max - panel_min,
    .groups = "drop"
  )

figure3_row_range_raw <- figure3_plot_data %>%
  mutate(row_panel = as.character(row_panel)) %>%
  group_by(row_panel) %>%
  summarise(
    row_min_raw = min(c(lwr, mean, 0), na.rm = TRUE),
    row_max_raw = max(c(upr, mean, 0), na.rm = TRUE),
    row_span_raw = row_max_raw - row_min_raw,
    .groups = "drop"
  )

# Prepare matching probabilities for the three capacity rows.
figure3_pmatch_base <- bind_rows(
  req_period_summary %>%
    transmute(
      scenario_lab = as.character(scenario_lab),
      period_mid,
      period_end,
      row_panel = "Plasticity",
      target_component = "Plasticity",
      Pmatch = mean_P_plasticity_match
    ),
  req_period_summary %>%
    transmute(
      scenario_lab = as.character(scenario_lab),
      period_mid,
      period_end,
      row_panel = "Evolution",
      target_component = "Evolution",
      Pmatch = mean_P_evolution_match
    ),
  req_period_summary %>%
    transmute(
      scenario_lab = as.character(scenario_lab),
      period_mid,
      period_end,
      row_panel = "Total capacity",
      target_component = "Total capacity",
      Pmatch = mean_P_match
    )
)

# Position labels above the highest error bar in each comparison.
figure3_pmatch_anchor <- figure3_pmatch_base %>%
  left_join(
    figure3_plot_data %>%
      mutate(
        scenario_lab = as.character(scenario_lab),
        row_panel = as.character(row_panel),
        component_lab = as.character(component_lab)
      ) %>%
      select(scenario_lab, period_mid, period_end, row_panel, component_lab, year_plot, mean, upr),
    by = c("scenario_lab", "period_mid", "period_end", "row_panel")
  ) %>%
  filter(component_lab == "Required" | component_lab == target_component) %>%
  group_by(scenario_lab, period_mid, period_end, row_panel, target_component, Pmatch) %>%
  summarise(
    anchor_x = mean(year_plot, na.rm = TRUE),
    anchor_y = max(mean, na.rm = TRUE),
    anchor_bar_top = max(upr, na.rm = TRUE),
    .groups = "drop"
  )

figure3_pmatch_data <- figure3_pmatch_anchor %>%
  left_join(
    figure3_panel_range,
    by = c("scenario_lab", "row_panel")
  ) %>%
  left_join(
    figure3_row_range_raw,
    by = "row_panel"
  ) %>%
  mutate(
    row_panel = factor(row_panel, levels = fig3_panel_levels),
    scenario_lab = factor(scenario_lab, levels = fig3_scenario_levels),
    label = sprintf("%.2f", Pmatch),
    label_x = anchor_x,
    label_y = anchor_bar_top + label_y_nudge * panel_span,
    label_hjust = 0.5,
    label_vjust = 0
  )

fig3_axis_theme <- theme_classic(base_size = 6.8) +
  theme(
    text = element_text(color = "#222222"),
    axis.title = element_blank(),
    axis.text = element_text(size = 7.2, color = "#222222"),
    axis.line = element_line(linewidth = 0.52, color = "#222222"),
    axis.ticks = element_line(linewidth = 0.46, color = "#222222"),
    axis.ticks.length = grid::unit(1.15, "mm"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.2, color = "#222222"),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    panel.spacing.x = grid::unit(0.55, "lines"),
    panel.spacing.y = grid::unit(0.55, "lines"),
    plot.margin = margin(1.8, 2.2, 1.8, 1.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

fig3_legend_labels <- c(
  "Required" = "required",
  "Plasticity" = "plasticity",
  "Evolution" = "evolution",
  "Total capacity" = "total"
)

# Use a common y-axis range within each row, including label positions.
figure3_row_limits <- bind_rows(
  figure3_plot_data %>%
    transmute(row_panel = as.character(row_panel), y_value = lwr),
  figure3_plot_data %>%
    transmute(row_panel = as.character(row_panel), y_value = upr),
  figure3_plot_data %>%
    transmute(row_panel = as.character(row_panel), y_value = mean),
  figure3_pmatch_data %>%
    transmute(row_panel = as.character(row_panel), y_value = label_y),
  data.frame(row_panel = fig3_panel_levels, y_value = 0)
) %>%
  group_by(row_panel) %>%
  summarise(
    y_min = min(y_value, na.rm = TRUE),
    y_max = max(y_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_span = y_max - y_min,
    y_min = y_min - 0.03 * y_span,
    y_max = y_max + 0.075 * y_span
  )

fig3_panel_theme <- fig3_axis_theme +
  theme(
    legend.position = "none",
    plot.margin = margin(1.2, 1.4, 1.2, 1.2)
  )

# Create one panel for each capacity and scenario.
make_figure3_panel <- function(
    row_value,
    scenario_value,
    show_x_axis = TRUE,
    show_y_axis = TRUE
) {
  panel_data <- figure3_plot_data %>%
    filter(
      as.character(row_panel) == row_value,
      as.character(scenario_lab) == scenario_value
    )

  panel_pmatch <- figure3_pmatch_data %>%
    filter(
      as.character(row_panel) == row_value,
      as.character(scenario_lab) == scenario_value
    )

  panel_limits <- figure3_row_limits %>%
    filter(row_panel == row_value)

  panel_plot <- ggplot(
    panel_data,
    aes(
      x = year_plot,
      y = mean,
      color = component_lab,
      group = component_lab
    )
  ) +
    geom_blank(aes(y = 0)) +
    geom_errorbar(
      aes(ymin = lwr, ymax = upr),
      width = 1.45,
      linewidth = 0.48,
      alpha = 0.78
    ) +
    geom_line(
      aes(linewidth = component_lab)
    ) +
    geom_point(
      size = 1.42,
      stroke = 0.22
    ) +
    geom_text(
      data = panel_pmatch,
      aes(
        x = label_x,
        y = label_y,
        label = label,
        hjust = label_hjust,
        vjust = label_vjust
      ),
      inherit.aes = FALSE,
      size = fig3_pmatch_label_size,
      fontface = "plain",
      color = "#222222"
    ) +
    scale_x_continuous(
      breaks = period_breaks,
      labels = period_labels_short,
      limits = fig3_x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 3),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      ylim = c(panel_limits$y_min, panel_limits$y_max),
      clip = "off"
    ) +
    scale_color_manual(
      values = fig3_cols,
      labels = fig3_legend_labels,
      drop = FALSE
    ) +
    scale_linewidth_manual(
      values = c(
        "Required" = 0.48,
        "Plasticity" = 0.62,
        "Evolution" = 0.62,
        "Total capacity" = 0.62
      ),
      guide = "none"
    ) +
    labs(x = NULL, y = NULL) +
    fig3_panel_theme

  if (!show_x_axis) {
    panel_plot <- panel_plot +
      theme(
        axis.text.x = element_text(size = 7.2, color = "transparent")
      )
  }

  if (!show_y_axis) {
    panel_plot <- panel_plot +
      theme(
        axis.text.y = element_text(size = 7.2, color = "transparent")
      )
  }

  panel_plot
}

# Show y-axis numbers in the first column and x-axis numbers in the last row.
figure3_panel_index <- bind_rows(lapply(seq_along(fig3_panel_levels), function(row_i) {
  data.frame(
    row_panel = fig3_panel_levels[row_i],
    scenario_lab = fig3_scenario_levels,
    row_index = row_i,
    col_index = seq_along(fig3_scenario_levels),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    panel_id = row_number(),
    show_y_axis = col_index == 1,
    show_x_axis = row_index == length(fig3_panel_levels)
  )

figure3_single_plots <- lapply(seq_len(nrow(figure3_panel_index)), function(i) {
  make_figure3_panel(
    row_value = figure3_panel_index$row_panel[i],
    scenario_value = figure3_panel_index$scenario_lab[i],
    show_x_axis = figure3_panel_index$show_x_axis[i],
    show_y_axis = figure3_panel_index$show_y_axis[i]
  )
})

# Create a shared legend.
figure3_legend_data <- data.frame(
  component_lab = factor(
    rep(names(fig3_cols), each = 2),
    levels = names(fig3_cols)
  ),
  x = rep(c(1, 2), times = length(fig3_cols)),
  y = rep(1, length(fig3_cols) * 2)
)

figure3_legend_plot <- ggplot(
  figure3_legend_data,
  aes(x = x, y = y, color = component_lab, group = component_lab)
) +
  geom_line(linewidth = 0.56) +
  geom_point(size = 1.25) +
  scale_color_manual(
    values = fig3_cols,
    labels = fig3_legend_labels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(linewidth = 0.56, size = 1.25)
    )
  ) +
  theme_void(base_size = 5.8) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.4, color = "#222222"),
    legend.key.height = grid::unit(2.3, "mm"),
    legend.key.width = grid::unit(3.3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

figure3_legend <- cowplot::get_legend(figure3_legend_plot)

# Set panel and legend export sizes.
fig3_panel_width <- 1.60
fig3_panel_height <- 1.18
fig3_legend_height <- 0.29
fig3_dpi <- if (exists("panel_size")) panel_size$dpi else 600

fig3_single_dir <- file.path(figure_dir, "figure3_panels")

if (!dir.exists(fig3_single_dir)) {
  dir.create(fig3_single_dir, recursive = TRUE)
}

clean_file_part <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  tolower(gsub("^_|_$", "", x))
}

# Save the 12 panels and the shared legend.
for (i in seq_len(nrow(figure3_panel_index))) {
  panel_stem <- paste0(
    sprintf("fig3_panel_%02d_", figure3_panel_index$panel_id[i]),
    clean_file_part(figure3_panel_index$row_panel[i]),
    "_",
    clean_file_part(figure3_panel_index$scenario_lab[i])
  )

  save_pub_pair(
    figure3_single_plots[[i]],
    file.path(fig3_single_dir, panel_stem),
    width = fig3_panel_width,
    height = fig3_panel_height,
    dpi = fig3_dpi
  )
}

save_pub_pair(
  cowplot::plot_grid(figure3_legend),
  file.path(fig3_single_dir, "fig3_shared_legend"),
  width = fig3_panel_width * 4,
  height = fig3_legend_height,
  dpi = fig3_dpi
)

# Assemble the three-by-four figure.
figure3_body <- cowplot::plot_grid(
  plotlist = figure3_single_plots,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

p_fig3 <- cowplot::plot_grid(
  figure3_body,
  figure3_legend,
  ncol = 1,
  rel_heights = c(1, 0.048)
)

print(p_fig3)

ggsave(
  file.path(figure_dir, "figure3.pdf"),
  plot = p_fig3,
  width = fig3_panel_width * 4,
  height = fig3_panel_height * 3 + 0.40,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(figure_dir, "figure3.png"),
  plot = p_fig3,
  width = fig3_panel_width * 4,
  height = fig3_panel_height * 3 + 0.40,
  units = "in",
  dpi = fig3_dpi,
  bg = "white"
)
