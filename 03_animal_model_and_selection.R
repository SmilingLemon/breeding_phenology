# Part 3: Random-regression animal model, selection, evolutionary response, and Figure 2
# Required packages: MCMCglmm, coda, dplyr, tidyr, ggplot2, scales.

public_dir <- "."
data_dir <- file.path(public_dir, "data")
output_dir <- file.path(public_dir, "results")
figure_dir <- file.path(output_dir, "figures")
data_file <- file.path(data_dir, "dxfamily.csv")
ped_file <- file.path(data_dir, "ped.csv")

# -----------------------------------------------------------------------------
# Part 0: Load packages and create output folders.
# -----------------------------------------------------------------------------

library(MCMCglmm)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Part 1: Prepare breeding data and the pedigree relationship matrix.
# -----------------------------------------------------------------------------

dx_data_rn <- read.csv(
  data_file,
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ped <- read.csv(
  ped_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# Prepare pedigree IDs and calculate the inverse relationship matrix.
colnames(ped) <- c("animal", "dam", "sire")

ped$animal <- as.character(ped$animal)
ped$dam <- as.character(ped$dam)
ped$sire <- as.character(ped$sire)

invAmatrix <- inverseA(ped)$Ainv

# Set numeric predictors and grouping factors.
dx_data_rn$ld41 <- as.numeric(dx_data_rn$ld41)
dx_data_rn$tem11_1 <- as.numeric(dx_data_rn$tem11_1)
dx_data_rn$age <- as.numeric(dx_data_rn$age)
dx_data_rn$TQ <- as.numeric(dx_data_rn$TQ)

dx_data_rn$animal <- as.factor(dx_data_rn$animal)
dx_data_rn$PI <- as.factor(dx_data_rn$PI)
dx_data_rn$year <- as.factor(dx_data_rn$year)

# Mean-center temperature, relative age, and territory quality without rescaling; the _z variables are centered, not standardized.
dx_data_rn$tem11_1_z <- as.numeric(scale(dx_data_rn$tem11_1, center = TRUE, scale = FALSE))
dx_data_rn$age_z <- as.numeric(scale(dx_data_rn$age, center = TRUE, scale = FALSE))
dx_data_rn$TQ_c <- as.numeric(scale(dx_data_rn$TQ, center = TRUE, scale = FALSE))

cat("========== Input checks ==========\n")
cat("Breeding records: ", nrow(dx_data_rn), "\n")
cat("Females: ", length(unique(dx_data_rn$female)), "\n")
cat("Animal IDs: ", length(unique(dx_data_rn$animal)), "\n")
cat("Permanent-environment IDs: ", length(unique(dx_data_rn$PI)), "\n")
cat("Years: ", length(unique(dx_data_rn$year)), "\n")
cat("Pedigree individuals: ", nrow(ped), "\n")
cat("mean tem11_1_z =", round(mean(dx_data_rn$tem11_1_z, na.rm = TRUE), 8), "\n")
cat("mean age_z =", round(mean(dx_data_rn$age_z, na.rm = TRUE), 8), "\n")
cat("mean TQ_c =", round(mean(dx_data_rn$TQ_c, na.rm = TRUE), 8), "\n\n")

# -----------------------------------------------------------------------------
# Part 2: Fit the random-regression animal model.
# -----------------------------------------------------------------------------

# Use 3,500,000 iterations, 500,000 burn-in, and thinning of 1,000 (3,000 retained samples).
Thin <- 1000
Burnin <- 500000
Nitt <- 3000 * Thin + Burnin

cat("========== MCMC settings ==========\n")
cat("Thin =", Thin, "\n")
cat("Burnin =", Burnin, "\n")
cat("Nitt =", Nitt, "\n\n")

# Animal model with additive-genetic, permanent-environment, and year effects.
prior_GE <- list(
  R = list(V = 1, nu = 1),
  G = list(
    G1 = list(
      V = diag(2),
      nu = 2,
      alpha.mu = rep(0, 2),
      alpha.V = diag(1000, 2, 2)
    ),
    G2 = list(
      V = diag(2),
      nu = 2,
      alpha.mu = rep(0, 2),
      alpha.V = diag(1000, 2, 2)
    ),
    G3 = list(
      V = 1,
      nu = 1,
      alpha.mu = 0,
      alpha.V = 1000
    )
  )
)

mod_GE <- MCMCglmm(
  fixed = ld41 ~ TQ_c + tem11_1_z + age_z,
  random = ~ us(1 + tem11_1_z):animal + us(1 + tem11_1_z):PI + year,
  rcov = ~ units,
  prior = prior_GE,
  nitt = Nitt,
  burnin = Burnin,
  thin = Thin,
  ginverse = list(animal = invAmatrix),
  data = dx_data_rn,
  pr = TRUE,
  saveX = TRUE,
  saveZ = TRUE,
  verbose = FALSE
)

summary(mod_GE)

saveRDS(
  mod_GE,
  file = file.path(output_dir, "mod_GE.rds")
)
cat("\nRandom-regression animal model saved.\n")

# -----------------------------------------------------------------------------
# Part 3: Summarize fixed effects.
# -----------------------------------------------------------------------------

# Read the saved model to summarize results without refitting.
mod_GE <- readRDS(file.path(output_dir, "mod_GE.rds"))

require_columns <- function(available, required, object_label) {
  missing_cols <- setdiff(required, available)

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        object_label, " is missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

fixed_GE <- summary(mod_GE)$solutions

# Order fixed effects as intercept, temperature, relative age, and territory quality.
required_fixed_terms <- c("(Intercept)", "tem11_1_z", "age_z", "TQ_c")

require_columns(
  rownames(fixed_GE),
  required_fixed_terms,
  "mod_GE fixed effects"
)

fixed_GE_reordered <- fixed_GE[required_fixed_terms, ]

fixed_table_GE <- data.frame(
  Model = "GxE with pedigree",
  Term = c("Intercept", "Winter temperature", "Female age", "Territory quality"),
  Estimate = fixed_GE_reordered[, "post.mean"],
  HPD_lower = fixed_GE_reordered[, "l-95% CI"],
  HPD_upper = fixed_GE_reordered[, "u-95% CI"],
  eff.samp = fixed_GE_reordered[, "eff.samp"],
  pMCMC = fixed_GE_reordered[, "pMCMC"],
  row.names = NULL
)

print(fixed_table_GE)

write.csv(
  fixed_table_GE,
  file.path(output_dir, "part3_mod_GE_fixed_effects.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 4: Summarize variance components and intercept heritability.
# -----------------------------------------------------------------------------

# Summarize posterior distributions with 95% HPD intervals.
get_hpd <- function(x) {
  HPDinterval(as.mcmc(x), prob = 0.95)
}

make_posterior_summary <- function(x, parameter, model) {
  h <- get_hpd(x)

  data.frame(
    Model = model,
    Parameter = parameter,
    Post.mean = mean(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    HPD_lower = as.numeric(h[1]),
    HPD_upper = as.numeric(h[2]),
    eff.samp = as.numeric(effectiveSize(as.mcmc(x))),
    row.names = NULL
  )
}

required_GE_vcv <- c(
  "(Intercept):(Intercept).animal",
  "tem11_1_z:tem11_1_z.animal",
  "tem11_1_z:(Intercept).animal",
  "(Intercept):(Intercept).PI",
  "tem11_1_z:tem11_1_z.PI",
  "tem11_1_z:(Intercept).PI",
  "year",
  "units"
)

require_columns(colnames(mod_GE$VCV), required_GE_vcv, "mod_GE$VCV")

VA_int_samples <- as.numeric(mod_GE$VCV[, "(Intercept):(Intercept).animal"])
VA_slope_samples <- as.numeric(mod_GE$VCV[, "tem11_1_z:tem11_1_z.animal"])
CoVA_samples <- as.numeric(mod_GE$VCV[, "tem11_1_z:(Intercept).animal"])

VPE_int_samples <- as.numeric(mod_GE$VCV[, "(Intercept):(Intercept).PI"])
VPE_slope_samples <- as.numeric(mod_GE$VCV[, "tem11_1_z:tem11_1_z.PI"])
CoVPE_samples <- as.numeric(mod_GE$VCV[, "tem11_1_z:(Intercept).PI"])

Vyear_samples <- as.numeric(mod_GE$VCV[, "year"])
Vresidual_samples <- as.numeric(mod_GE$VCV[, "units"])

# Calculate intercept heritability at the centered temperature.
VP_int_samples <- VA_int_samples +
  VPE_int_samples +
  Vyear_samples +
  Vresidual_samples

h2_int_samples <- VA_int_samples / VP_int_samples
h2_int_HPD <- get_hpd(h2_int_samples)

variance_summary <- bind_rows(
  make_posterior_summary(VA_int_samples, "Additive genetic variance, intercept", "GxE with pedigree"),
  make_posterior_summary(VA_slope_samples, "Additive genetic variance, slope", "GxE with pedigree"),
  make_posterior_summary(CoVA_samples, "Additive genetic covariance", "GxE with pedigree"),
  make_posterior_summary(VPE_int_samples, "Permanent environmental variance, intercept", "GxE with pedigree"),
  make_posterior_summary(VPE_slope_samples, "Permanent environmental variance, slope", "GxE with pedigree"),
  make_posterior_summary(CoVPE_samples, "Permanent environmental covariance", "GxE with pedigree"),
  make_posterior_summary(Vyear_samples, "Year variance", "GxE with pedigree"),
  make_posterior_summary(Vresidual_samples, "Residual variance", "GxE with pedigree"),
  make_posterior_summary(h2_int_samples, "Intercept heritability", "GxE with pedigree")
)

print(variance_summary)

write.csv(
  variance_summary,
  file.path(output_dir, "part3_variance_component_summary.csv"),
  row.names = FALSE
)

random_table_GE <- variance_summary %>%
  filter(Model == "GxE with pedigree") %>%
  mutate(
    Estimate = round(Post.mean, 2),
    HPD = paste0(round(HPD_lower, 2), ", ", round(HPD_upper, 2)),
    eff.samp = round(eff.samp, 0)
  ) %>%
  select(Parameter, Estimate, HPD, eff.samp)
print(random_table_GE)

write.csv(
  random_table_GE,
  file.path(output_dir, "part3_mod_GE_random_effects.csv"),
  row.names = FALSE
)

cat("h2_int mean =", round(mean(h2_int_samples), 4), "\n")
cat("h2_int 95% HPD =", round(h2_int_HPD[1], 4), ",", round(h2_int_HPD[2], 4), "\n")

# -----------------------------------------------------------------------------
# Part 5: Selection on raw reaction norms and evolutionary response.
# -----------------------------------------------------------------------------

if (!exists("h2_int_samples")) {
  stop("h2_int_samples was not found. Please run Part 4 first.")
}

# Use females with at least two breeding records for raw reaction norms.
if (exists("dx_data_rn")) {
  dx_data_sel <- dx_data_rn
} else {
  if (!exists("data_file")) {
    stop("data_file was not found. Please run Part 1 first.")
  }

  dx_data_sel <- read.csv(
    data_file,
    fileEncoding = "UTF-8-BOM",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

required_selection_cols <- c("female", "fledgling", "ld41", "tem11_1")
require_columns(names(dx_data_sel), required_selection_cols, "selection data")

fledgling_char <- as.character(dx_data_sel$fledgling)
fledgling_char[fledgling_char == "#N/A"] <- NA
dx_data_sel$fledgling <- as.numeric(fledgling_char)

dx_data_sel$ld41 <- as.numeric(dx_data_sel$ld41)
dx_data_sel$tem11_1 <- as.numeric(dx_data_sel$tem11_1)

if (!"tem11_1_z" %in% names(dx_data_sel)) {
  dx_data_sel$tem11_1_z <- as.numeric(
    scale(dx_data_sel$tem11_1, center = TRUE, scale = FALSE)
  )
}

raw_dat <- dx_data_sel %>%
  filter(
    !is.na(female),
    !is.na(ld41),
    !is.na(tem11_1_z)
  )

female_keep <- raw_dat %>%
  group_by(female) %>%
  summarise(
    n_record = n(),
    .groups = "drop"
  ) %>%
  filter(n_record >= 2)

raw_dat_keep <- raw_dat %>%
  filter(female %in% female_keep$female)

# Fit a separate laying-date regression on centered temperature for each repeat female.
raw_rn_df <- raw_dat_keep %>%
  group_by(female) %>%
  group_modify(~ {
    fit <- lm(ld41 ~ tem11_1_z, data = .x)

    data.frame(
      raw_intercept = coef(fit)[["(Intercept)"]],
      raw_slope = coef(fit)[["tem11_1_z"]]
    )
  }) %>%
  ungroup()

# Scale each female mean fledgling number by the mean across all 302 females before selecting repeat females.
fitness_data <- dx_data_sel %>%
  group_by(female) %>%
  summarise(
    fledgling = mean(fledgling, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rel_fitness = fledgling / mean(fledgling, na.rm = TRUE)
  )

df_sel <- raw_rn_df %>%
  inner_join(fitness_data, by = "female") %>%
  filter(
    !is.na(raw_intercept),
    !is.na(raw_slope),
    !is.na(rel_fitness)
  )

if (nrow(df_sel) < 3) {
  stop("Too few females are available for the selection analysis.")
}

cat("Total female count:", length(unique(dx_data_sel$female)), "\n")
cat("Females with >= 2 breeding records:", nrow(female_keep), "\n")
cat("Females in selection analysis:", nrow(df_sel), "\n")
cat("mean raw_intercept =", round(mean(df_sel$raw_intercept), 3), "\n")
cat("mean raw_slope =", round(mean(df_sel$raw_slope), 3), "\n")

write.csv(
  df_sel,
  file.path(output_dir, "part3_selection_reaction_norm_data.csv"),
  row.names = FALSE
)

# Fit unweighted fitness regressions on the raw intercepts and slopes.
selection_model_intercept <- lm(rel_fitness ~ raw_intercept, data = df_sel)
selection_model_slope <- lm(rel_fitness ~ raw_slope, data = df_sel)

get_lm_result <- function(model, term, label) {
  s <- summary(model)$coefficients
  ci <- confint(model)

  data.frame(
    Trait = label,
    Term = term,
    beta = s[term, "Estimate"],
    SE = s[term, "Std. Error"],
    lower_95CI = ci[term, 1],
    upper_95CI = ci[term, 2],
    P_value = s[term, "Pr(>|t|)"],
    row.names = NULL
  )
}

selection_results <- bind_rows(
  get_lm_result(
    selection_model_intercept,
    "raw_intercept",
    "Reaction-norm intercept"
  ),
  get_lm_result(
    selection_model_slope,
    "raw_slope",
    "Reaction-norm slope"
  )
)

print(selection_results)

write.csv(
  selection_results,
  file.path(output_dir, "part3_selection_regression_results.csv"),
  row.names = FALSE
)

# Bootstrap the intercept selection differential using 5,000 female-level resamples.
set.seed(123)

S_est <- cov(df_sel$raw_intercept, df_sel$rel_fitness, use = "complete.obs")

B_S <- 5000
valid_rows <- which(complete.cases(df_sel$raw_intercept, df_sel$rel_fitness))

if (length(valid_rows) < 2) {
  stop("Too few complete rows are available to bootstrap S.")
}

S_boot <- replicate(B_S, {
  rows_sel <- sample(valid_rows, length(valid_rows), replace = TRUE)
  cov(df_sel$raw_intercept[rows_sel], df_sel$rel_fitness[rows_sel])
})

S_CI <- quantile(S_boot, probs = c(0.025, 0.975), na.rm = TRUE)

cat("S =", round(S_est, 4), "\n")
cat("S 95% bootstrap CI =", round(S_CI[1], 4), ",", round(S_CI[2], 4), "\n")

# Propagate bootstrap selection and posterior heritability draws through R = h2 * S.
set.seed(123)

n_draw <- min(length(h2_int_samples), length(S_boot))

cat("h2_int_samples length =", length(h2_int_samples), "\n")
cat("S_boot length =", length(S_boot), "\n")
cat("n_draw used for R =", n_draw, "\n")

h2_sample <- sample(h2_int_samples, n_draw, replace = TRUE)
S_sample <- sample(S_boot, n_draw, replace = TRUE)

R_post <- h2_sample * S_sample
R_HPD <- HPDinterval(as.mcmc(R_post), prob = 0.95)
h2_int_HPD <- HPDinterval(as.mcmc(h2_int_samples), prob = 0.95)

cat("R_post mean =", round(mean(R_post), 4), "days/generation\n")
cat("R_post 95% HPD =", round(R_HPD[1], 4), ",", round(R_HPD[2], 4), "days/generation\n")

ER_table <- data.frame(
  Parameter = c(
    "Selection differential (S)",
    "Intercept heritability (h2_int)",
    "Per-generation evolutionary response (R)"
  ),
  Estimate = as.numeric(c(
    S_est,
    mean(h2_int_samples),
    mean(R_post)
  )),
  lower_95 = as.numeric(c(
    S_CI[1],
    h2_int_HPD[1],
    R_HPD[1]
  )),
  upper_95 = as.numeric(c(
    S_CI[2],
    h2_int_HPD[2],
    R_HPD[2]
  )),
  Interval_type = c(
    "Bootstrap 95% CI",
    "Bayesian 95% HPD",
    "Propagated 95% uncertainty interval"
  ),
  n_samples = c(
    length(S_boot),
    length(h2_int_samples),
    length(R_post)
  ),
  Effective_sample_size = c(
    NA,
    as.numeric(effectiveSize(as.mcmc(h2_int_samples))),
    as.numeric(effectiveSize(as.mcmc(R_post)))
  )
)

print(ER_table)

write.csv(
  ER_table,
  file.path(output_dir, "part3_evolutionary_response_summary.csv"),
  row.names = FALSE
)

evolutionary_response_samples <- data.frame(
  h2_int = h2_sample,
  S = S_sample,
  R = R_post
)

write.csv(
  evolutionary_response_samples,
  file.path(output_dir, "part3_evolutionary_response_posterior_samples.csv"),
  row.names = FALSE
)
saveRDS(
  evolutionary_response_samples,
  file.path(output_dir, "part3_evolutionary_response_posterior_samples.rds")
)

# -----------------------------------------------------------------------------
# Part 6: Figures 2a-b, posterior mean reaction norms.
# -----------------------------------------------------------------------------

source(file.path(public_dir, "plot_style_publication.R"))

# Set reaction-norm colors.
fig2_palette <- c(
  main_blue = "#4E79A7",
  ribbon_blue = "#D7E3EE",
  muted_purple = "#8E77B8",
  muted_teal = "#59A89C",
  muted_orange = "#D99045",
  phen_blue = "#4CA8BF",
  black = "#222222",
  fig2_genetic_color <- "#F6949A"
)

# Read breeding records for the phenotyped females.
dx_data_rn <- read.csv(
  data_file,
  fileEncoding = "UTF-8-BOM",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dx_data_rn <- dx_data_rn %>%
  mutate(
    animal = as.character(animal),
    PI = as.character(PI),
    female = as.character(female),
    ld41 = as.numeric(ld41),
    tem11_1 = as.numeric(tem11_1),
    age = as.numeric(age),
    fledgling = na_if(as.character(fledgling), "#N/A"),
    fledgling = as.numeric(fledgling)
  )

if (!"tem11_1_z" %in% names(dx_data_rn)) {
  dx_data_rn$tem11_1_z <- as.numeric(
    scale(dx_data_rn$tem11_1, center = TRUE, scale = FALSE)
  )
}

cat("Breeding records:", nrow(dx_data_rn), "\n")
cat("mod_GE posterior samples:", nrow(mod_GE$Sol), "\n\n")

# Extract posterior mean animal and permanent-environment effects.
sol_i <- colMeans(mod_GE$Sol)

animal_int <- sol_i[grep("^\\(Intercept\\)\\.animal\\.", names(sol_i))]
animal_slope <- sol_i[grep("^tem11_1_z\\.animal\\.", names(sol_i))]
pi_int <- sol_i[grep("^\\(Intercept\\)\\.PI\\.", names(sol_i))]
pi_slope <- sol_i[grep("^tem11_1_z\\.PI\\.", names(sol_i))]

animal_ids <- sub("^\\(Intercept\\)\\.animal\\.", "", names(animal_int))
animal_slope_ids <- sub("^tem11_1_z\\.animal\\.", "", names(animal_slope))
pi_int_ids <- sub("^\\(Intercept\\)\\.PI\\.", "", names(pi_int))
pi_slope_ids <- sub("^tem11_1_z\\.PI\\.", "", names(pi_slope))

reaction_effects <- data.frame(
  id = animal_ids,
  a_int = as.numeric(animal_int),
  a_slope = as.numeric(animal_slope[match(animal_ids, animal_slope_ids)]),
  pi_int = as.numeric(pi_int[match(animal_ids, pi_int_ids)]),
  pi_slope = as.numeric(pi_slope[match(animal_ids, pi_slope_ids)]),
  stringsAsFactors = FALSE
)

phenotyped_animal <- unique(dx_data_rn$animal[!is.na(dx_data_rn$animal)])

reaction_effects <- reaction_effects %>%
  filter(id %in% phenotyped_animal)

# Evaluate deviations across centered temperatures from -2 to 2 degrees Celsius.
E_seq <- seq(-2, 2, length.out = 100)

reaction_plot_data <- expand.grid(
  id = reaction_effects$id,
  tem11_1_z = E_seq
) %>%
  left_join(reaction_effects, by = "id") %>%
  mutate(
# Genetic deviation: animal intercept and slope.
    y_gen = a_int + a_slope * tem11_1_z,
# Phenotypic deviation: animal plus permanent-environment intercepts and slopes.
    y_phen = (a_int + pi_int) + (a_slope + pi_slope) * tem11_1_z
  )

reaction_plot_phen <- reaction_plot_data %>%
  filter(!is.na(y_phen))

reaction_plot_gen <- reaction_plot_data %>%
  filter(!is.na(y_gen))

# Use the same y-axis range for both reaction-norm panels.
y_lim_reaction <- range(
  c(reaction_plot_data$y_phen, reaction_plot_data$y_gen),
  na.rm = TRUE
)

# Figure 2a: animal plus permanent-environment effects.
p_fig2a <- ggplot(
  reaction_plot_phen,
  aes(x = tem11_1_z, y = y_phen, group = id)
) +
  geom_line(
    color = fig2_palette[["phen_blue"]],
    alpha = 0.30,
    linewidth = 0.24
  ) +
  coord_cartesian(ylim = y_lim_reaction, clip = "off") +
  scale_x_continuous(breaks = seq(-2, 2, by = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8)

# Figure 2b: animal effects.
p_fig2b <- ggplot(
  reaction_plot_gen,
  aes(x = tem11_1_z, y = y_gen, group = id)
) +
  geom_line(
    color = fig2_genetic_color,
    alpha = 0.34,
    linewidth = 0.25
  ) +
  coord_cartesian(ylim = y_lim_reaction, clip = "off") +
  scale_x_continuous(breaks = seq(-2, 2, by = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8)

print(p_fig2a)
print(p_fig2b)

write.csv(
  reaction_plot_phen,
  file.path(output_dir, "figure2a_phenotypic_reaction_norm_data.csv"),
  row.names = FALSE
)

write.csv(
  reaction_plot_gen,
  file.path(output_dir, "figure2b_genetic_reaction_norm_data.csv"),
  row.names = FALSE
)

save_pub_pair(
  fix_pub_panel_size(p_fig2a),
  file.path(figure_dir, "fig2a_phenotypic_reaction_norms"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

save_pub_pair(
  fix_pub_panel_size(p_fig2b),
  file.path(figure_dir, "fig2b_genetic_reaction_norms"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigures 2a-b saved to results/figures.\n")

# -----------------------------------------------------------------------------
# Part 7: Figure 2c, posterior densities of intercept variances.
# -----------------------------------------------------------------------------

# Compare additive-genetic and permanent-environment intercept variance draws.
fig2c_data <- data.frame(
  value = c(VA_int_samples, VPE_int_samples),
  component = factor(
    rep(
      c("Genetic", "Permanent environment"),
      each = length(VA_int_samples)
    ),
    levels = c("Genetic", "Permanent environment")
  )
)

density_cols <- c(
  "Genetic" = "#F6949A",
  "Permanent environment" = "#2FADAE"
)

p_fig2c <- ggplot(
  fig2c_data,
  aes(x = value, fill = component, color = component)
) +
  geom_density(
    alpha = 0.42,
    linewidth = 0.42,
    adjust = 1.5
  ) +
  scale_fill_manual(values = density_cols) +
  scale_color_manual(values = density_cols) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8) +
  theme(
    legend.position = "none",
    legend.justification = c(0.5, 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.0),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.width = grid::unit(3.4, "mm")
  )

print(p_fig2c)

save_pub_pair(
  fix_pub_panel_size(p_fig2c),
  file.path(figure_dir, "fig2c_intercept_variance_density"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

# -----------------------------------------------------------------------------
# Part 8: Figure 2d, posterior densities of slope variances.
# -----------------------------------------------------------------------------

# Compare additive-genetic and permanent-environment slope variance draws.
fig2d_data <- data.frame(
  value = c(VA_slope_samples, VPE_slope_samples),
  component = factor(
    rep(
      c("Genetic", "Permanent environment"),
      each = length(VA_slope_samples)
    ),
    levels = c("Genetic", "Permanent environment")
  )
)

p_fig2d <- ggplot(
  fig2d_data,
  aes(x = value, fill = component, color = component)
) +
  geom_density(
    alpha = 0.42,
    linewidth = 0.42,
    adjust = 1.5
  ) +
  scale_fill_manual(values = density_cols) +
  scale_color_manual(values = density_cols) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8) +
  theme(
    legend.position = "none"
  )

print(p_fig2d)

save_pub_pair(
  fix_pub_panel_size(p_fig2d),
  file.path(figure_dir, "fig2d_slope_variance_density"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

# -----------------------------------------------------------------------------
# Part 9: Figures 2e-f, posterior mean animal effects and relative fitness.
# -----------------------------------------------------------------------------

# Set point, line, and regression-ribbon styles.
fig2_scatter <- list(
  color = fig2_palette[["main_blue"]],
  ribbon = fig2_palette[["ribbon_blue"]],
  point_size = pub_scatter$point_size,
  point_stroke = pub_scatter$point_stroke,
  line_width = pub_scatter$line_width
)

if (!exists("dx_data_rn")) {
  dx_data_rn <- read.csv(
    data_file,
    fileEncoding = "UTF-8-BOM",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

dx_data_rn <- dx_data_rn %>%
  mutate(
    female = as.character(female),
    fledgling = na_if(as.character(fledgling), "#N/A"),
    fledgling = as.numeric(fledgling)
  )

# Extract posterior mean animal intercepts and slopes for the phenotyped females.
BV_cols <- grep("animal", colnames(mod_GE$Sol), value = TRUE)
BV_int_cols <- grep("\\(Intercept\\).*\\.animal", BV_cols, value = TRUE)
BV_slope_cols <- grep("tem11_1_z.*\\.animal", BV_cols, value = TRUE)

if (length(BV_int_cols) == 0 || length(BV_slope_cols) == 0) {
  stop("Could not find animal intercept or slope columns in mod_GE$Sol for Figure 2e-f.")
}

intercepts <- colMeans(mod_GE$Sol[, BV_int_cols, drop = FALSE])
slopes <- colMeans(mod_GE$Sol[, BV_slope_cols, drop = FALSE])

names(intercepts) <- gsub(".*\\.animal\\.", "", names(intercepts))
names(slopes) <- gsub(".*\\.animal\\.", "", names(slopes))

all_BV_df <- data.frame(
  female = names(intercepts),
  intercept = as.numeric(intercepts),
  slope = as.numeric(slopes[names(intercepts)]),
  stringsAsFactors = FALSE
)

reaction_norm_df <- all_BV_df %>%
  filter(female %in% dx_data_rn$female)

fitness_data <- dx_data_rn %>%
  group_by(female) %>%
  summarise(
    fledgling = mean(fledgling, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rel_fitness = fledgling / mean(fledgling, na.rm = TRUE)
  )

df_sel <- reaction_norm_df %>%
  inner_join(fitness_data, by = "female") %>%
  filter(
    !is.na(intercept),
    !is.na(slope),
    !is.na(rel_fitness)
  ) %>%
  mutate(
# Standardize posterior mean intercepts and slopes across the females in this analysis.
    std_intercept = as.numeric(scale(intercept)),
    std_slope = as.numeric(scale(slope))
  )

if (nrow(df_sel) < 3) {
  stop("Too few females are available for the Figure 2e-f selection plots.")
}

selection_model_intercept <- lm(
  rel_fitness ~ std_intercept,
  data = df_sel
)

selection_model_slope <- lm(
  rel_fitness ~ std_slope,
  data = df_sel
)

print(summary(selection_model_intercept))
print(confint(selection_model_intercept))

print(summary(selection_model_slope))
print(confint(selection_model_slope))

cat("n =", nrow(df_sel), "females\n")

# Calculate fitted values and 95% confidence intervals for the regression lines.
pred_e <- data.frame(
  std_intercept = seq(
    min(df_sel$std_intercept, na.rm = TRUE),
    max(df_sel$std_intercept, na.rm = TRUE),
    length.out = 100
  )
)

pred_e <- cbind(
  pred_e,
  as.data.frame(predict(
    selection_model_intercept,
    newdata = pred_e,
    interval = "confidence"
  ))
)

pred_f <- data.frame(
  std_slope = seq(
    min(df_sel$std_slope, na.rm = TRUE),
    max(df_sel$std_slope, na.rm = TRUE),
    length.out = 100
  )
)

pred_f <- cbind(
  pred_f,
  as.data.frame(predict(
    selection_model_slope,
    newdata = pred_f,
    interval = "confidence"
  ))
)

write.csv(
  df_sel,
  file.path(figure_dir, "figure2ef_posterior_animal_effect_fitness_data.csv"),
  row.names = FALSE
)

# Figure 2e: standardized animal intercept and relative fitness.
p_fig2e <- ggplot(df_sel, aes(x = std_intercept, y = rel_fitness)) +
  geom_point(
    shape = 21,
    size = fig2_scatter$point_size,
    stroke = fig2_scatter$point_stroke,
    color = fig2_scatter$color,
    fill = "white",
    alpha = 0.88
  ) +
  geom_ribbon(
    data = pred_e,
    aes(x = std_intercept, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = pub_alpha(fig2_scatter$ribbon, 0.8)
  ) +
  geom_line(
    data = pred_e,
    aes(x = std_intercept, y = fit),
    inherit.aes = FALSE,
    color = fig2_scatter$color,
    linewidth = fig2_scatter$line_width
  )  +scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(x = NULL, y = NULL) +
  pub_theme(base_size = 5.8)

# Figure 2f: standardized animal slope and relative fitness.
p_fig2f <- ggplot(df_sel, aes(x = std_slope, y = rel_fitness)) +
  geom_point(
    shape = 21,
    size = fig2_scatter$point_size,
    stroke = fig2_scatter$point_stroke,
    color = fig2_scatter$color,
    fill = "white",
    alpha = 0.88
  ) +
  geom_ribbon(
    data = pred_f,
    aes(x = std_slope, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = pub_alpha(fig2_scatter$ribbon, 0.8)
  ) +
  geom_line(
    data = pred_f,
    aes(x = std_slope, y = fit),
    inherit.aes = FALSE,
    color = fig2_scatter$color,
    linewidth = fig2_scatter$line_width,
    linetype = "dashed"
  ) +scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(x = NULL, y = NULL) +
  pub_theme(base_size = 5.8)

print(p_fig2e)
print(p_fig2f)

save_pub_pair(
  fix_pub_panel_size(p_fig2e),
  file.path(figure_dir, "fig2e_intercept_relative_fitness"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

save_pub_pair(
  fix_pub_panel_size(p_fig2f),
  file.path(figure_dir, "fig2f_slope_relative_fitness"),
  width = panel_size$single_width,
  height = panel_size$single_height
)
