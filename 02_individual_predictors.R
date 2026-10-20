# Part 2: Individual-level predictors of breeding time and Figure 1f
# Required packages: MCMCglmm, coda, dplyr, ggplot2, scales.

public_dir <- "."
data_dir <- file.path(public_dir, "data")
output_dir <- file.path(public_dir, "results")
figure_dir <- file.path(output_dir, "figures")

# -----------------------------------------------------------------------------
# Part 0: Load packages and create output folders.
# -----------------------------------------------------------------------------

library(MCMCglmm)
library(coda)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Part 1: Prepare the breeding data.
# -----------------------------------------------------------------------------

dx_data <- read.csv(
  file.path(data_dir, "dxfamily.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove the UTF-8 BOM from column names.
names(dx_data) <- sub("^\xef\xbb\xbf", "", names(dx_data), useBytes = TRUE)

# Check that all model variables are present.
required_vars <- c("year", "female", "ld41", "tem11_1", "age", "TQ", "coop")
missing_vars <- setdiff(required_vars, names(dx_data))

if (length(missing_vars) > 0) {
  stop(paste("Missing variables: ", paste(missing_vars, collapse = ", ")))
}

# Set numeric predictors and grouping factors.
dx_data$ld41 <- as.numeric(dx_data$ld41)
dx_data$tem11_1 <- as.numeric(dx_data$tem11_1)
dx_data$age <- as.numeric(dx_data$age)
dx_data$TQ <- as.numeric(dx_data$TQ)
dx_data$coop <- as.numeric(dx_data$coop)

dx_data$female <- as.factor(dx_data$female)
dx_data$year <- as.factor(dx_data$year)

# Mean-center temperature, relative age, and territory quality without rescaling.
dx_data$tem11_1_c <- dx_data$tem11_1 - mean(dx_data$tem11_1, na.rm = TRUE)
dx_data$age_c <- dx_data$age - mean(dx_data$age, na.rm = TRUE)
dx_data$TQ_c <- dx_data$TQ - mean(dx_data$TQ, na.rm = TRUE)

# Check the binary helper-presence coding.
coop_values <- sort(unique(dx_data$coop[!is.na(dx_data$coop)]))

if (!all(coop_values %in% c(0, 1))) {
  stop("coop must contain only 0 and 1, or missing values.")
}

# Retain complete records for the model variables.
model_vars <- c(
  "ld41",
  "tem11_1_c",
  "age_c",
  "TQ_c",
  "coop",
  "female",
  "year"
)

model_data <- dx_data[complete.cases(dx_data[, model_vars]), model_vars]
model_data <- droplevels(model_data)

cat("========== Input checks ==========\n")
cat("Input records: ", nrow(dx_data), "\n")
cat("Model records: ", nrow(model_data), "\n")
cat("Females: ", length(unique(model_data$female)), "\n")
cat("Years: ", length(unique(model_data$year)), "\n")
cat("coop values: ", paste(coop_values, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Part 2: Set variance priors and MCMC parameters.
# -----------------------------------------------------------------------------

# Variance priors for residuals and the female and year random intercepts.
prior_individual <- list(
  R = list(V = 1, nu = 1),
  G = list(
    G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
    G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)
  )
)

Thin <- 1000
Burnin <- 500000
Nitt <- 3000 * Thin + Burnin

cat("========== MCMC settings ==========\n")
cat("Thin =", Thin, "\n")
cat("Burnin =", Burnin, "\n")
cat("Nitt =", Nitt, "\n\n")

# -----------------------------------------------------------------------------
# Part 3: Fit the Bayesian Gaussian mixed model.
# -----------------------------------------------------------------------------

# Use 3,500,000 iterations, 500,000 burn-in, and thinning of 1,000 (3,000 retained samples).
set.seed(123)

part2_bayesian_individual_model <- MCMCglmm(
  fixed = ld41 ~ tem11_1_c + age_c + TQ_c + coop,
  random = ~ female + year,
  rcov = ~ units,
  family = "gaussian",
  prior = prior_individual,
  data = model_data,
  nitt = Nitt,
  burnin = Burnin,
  thin = Thin,
  verbose = FALSE
)

summary(part2_bayesian_individual_model)

saveRDS(
  part2_bayesian_individual_model,
  file = file.path(output_dir, "part2_bayesian_individual_model.rds")
)

cat("\nModel saved to results/part2_bayesian_individual_model.rds\n")

# -----------------------------------------------------------------------------
# Part 4: Summarize fixed effects with 95% HPD intervals.
# -----------------------------------------------------------------------------

# Read the saved model to summarize results without refitting.
part2_bayesian_individual_model <- readRDS(file.path(output_dir, "part2_bayesian_individual_model.rds"))

model_summary <- summary(part2_bayesian_individual_model)

# Summarize posterior means, standard deviations, 95% HPD intervals, and effective sample sizes.
sol_matrix <- as.matrix(part2_bayesian_individual_model$Sol)
sol_ess <- effectiveSize(part2_bayesian_individual_model$Sol)
sol_hpd <- HPDinterval(part2_bayesian_individual_model$Sol, prob = 0.95)

fixed_effects <- data.frame(
  Term = colnames(sol_matrix),
  Posterior_mean = colMeans(sol_matrix),
  Posterior_SD = apply(sol_matrix, 2, sd),
  HPD_lower = sol_hpd[colnames(sol_matrix), "lower"],
  HPD_upper = sol_hpd[colnames(sol_matrix), "upper"],
  Effective_sample_size = as.numeric(sol_ess[colnames(sol_matrix)]),
  pMCMC = model_summary$solutions[colnames(sol_matrix), "pMCMC"],
  row.names = NULL
)

fixed_effects$Term_label <- fixed_effects$Term
fixed_effects$Term_label[fixed_effects$Term == "(Intercept)"] <- "Intercept"
fixed_effects$Term_label[fixed_effects$Term == "tem11_1_c"] <- "Winter temperature"
fixed_effects$Term_label[fixed_effects$Term == "age_c"] <- "Age"
fixed_effects$Term_label[fixed_effects$Term == "TQ_c"] <- "Territory quality"
fixed_effects$Term_label[fixed_effects$Term == "coop"] <- "Cooperation"

fixed_effects <- fixed_effects[, c(
  "Term",
  "Term_label",
  "Posterior_mean",
  "Posterior_SD",
  "HPD_lower",
  "HPD_upper",
  "Effective_sample_size",
  "pMCMC"
)]

print(fixed_effects)
cat("\n")

write.csv(
  fixed_effects,
  file = file.path(output_dir, "part2_bayesian_individual_fixed_effects.csv"),
  row.names = FALSE
)

cat("Fixed effects saved to results/part2_bayesian_individual_fixed_effects.csv\n")

# -----------------------------------------------------------------------------
# Part 5: Summarize random-intercept and residual variances.
# -----------------------------------------------------------------------------

vcv_matrix <- as.matrix(part2_bayesian_individual_model$VCV)
vcv_ess <- effectiveSize(part2_bayesian_individual_model$VCV)
vcv_hpd <- HPDinterval(part2_bayesian_individual_model$VCV, prob = 0.95)

random_effects <- data.frame(
  Term = colnames(vcv_matrix),
  Posterior_mean = colMeans(vcv_matrix),
  Posterior_SD = apply(vcv_matrix, 2, sd),
  HPD_lower = vcv_hpd[colnames(vcv_matrix), "lower"],
  HPD_upper = vcv_hpd[colnames(vcv_matrix), "upper"],
  Effective_sample_size = as.numeric(vcv_ess[colnames(vcv_matrix)]),
  row.names = NULL
)

random_effects$Term_label <- random_effects$Term
random_effects$Term_label[random_effects$Term == "female"] <- "Female identity"
random_effects$Term_label[random_effects$Term == "year"] <- "Year"
random_effects$Term_label[random_effects$Term == "units"] <- "Residual"

random_effects <- random_effects[, c(
  "Term",
  "Term_label",
  "Posterior_mean",
  "Posterior_SD",
  "HPD_lower",
  "HPD_upper",
  "Effective_sample_size"
)]

print(random_effects)

write.csv(
  random_effects,
  file = file.path(output_dir, "part2_bayesian_individual_random_effects.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 6: Figure 1f, posterior means and 95% HPD intervals.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

source(file.path(public_dir, "plot_style_publication.R"))

fig1f_color <- "#4E79A7"

fixed_effects <- read.csv(
  file.path(output_dir, "part2_bayesian_individual_fixed_effects.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Plot the four predictors, excluding the intercept.
fig1f_data <- fixed_effects %>%
  filter(Term != "(Intercept)") %>%
  mutate(
    Term_label = case_when(
      Term == "tem11_1_c" ~ "Tmp",
      Term == "age_c" ~ "Age",
      Term == "TQ_c" ~ "TQ",
      Term == "coop" ~ "Co",
      TRUE ~ Term
    ),
    Term_label = factor(
      Term_label,
      levels = c("Co", "TQ", "Age", "Tmp")
    ),
    crosses_zero = HPD_lower <= 0 & HPD_upper >= 0
  )

fig1f_cols <- c(
  "Tmp" = fig1f_color,
  "Age" = fig1f_color,
  "TQ" = fig1f_color,
  "Co" = fig1f_color
)

fig1f_x_range <- range(
  c(fig1f_data$HPD_lower, fig1f_data$HPD_upper, fig1f_data$Posterior_mean),
  na.rm = TRUE
)

fig1f_x_span <- diff(fig1f_x_range)

# Set the x-axis limits from the displayed estimates and HPD intervals.
fig1f_x_limits <- c(
  fig1f_x_range[1] - 0.06 * fig1f_x_span,
  fig1f_x_range[2] + 0.06 * fig1f_x_span
)

cat("\n========== Figure 1f plotting data ==========\n")
print(fig1f_data)

# Points show posterior means; horizontal bars show 95% HPD intervals.
p_fig1f <- ggplot(
  fig1f_data,
  aes(y = Term_label, x = Posterior_mean, color = Term_label)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.25,
    color = pub_palette[["grey"]]
  ) +
  geom_errorbar(
    aes(xmin = HPD_lower, xmax = HPD_upper, alpha = crosses_zero),
    orientation = "y",
    width = 0.13,
    linewidth = 0.42
  ) +
  geom_point(
    aes(alpha = crosses_zero),
    size = 1.15
  ) +
  scale_color_manual(values = fig1f_cols) +
  scale_alpha_manual(
    values = c("TRUE" = 0.72, "FALSE" = 1.0),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 4),
    limits = fig1f_x_limits
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  pub_theme(base_size = 5.8) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.36, color = pub_palette[["black"]]),
    plot.margin = margin(2.5, 3.0, 2.5, 2.5)
  )

print(p_fig1f)

save_pub_pair(
  fix_pub_panel_size(p_fig1f),
  file.path(figure_dir, "fig1f_bayesian_individual_predictors"),
  width = panel_size$single_width,
  height = panel_size$single_height
)

cat("\nFigure 1f saved to results/figures.\n")
