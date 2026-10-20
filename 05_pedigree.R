# Table S1 pedigree summaries and annual climate trends.
# Required packages: dplyr, MCMCglmm, terra.
# Download WorldClim historical monthly tmin, tmax and prec (2.5 arc-minutes) separately.

library(dplyr)
library(MCMCglmm)
# Part 0: Set input and output paths.
public_dir <- "."
data_dir <- file.path(public_dir, "data")
output_dir <- file.path(public_dir, "results")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ped_file <- file.path(data_dir, "ped.csv")
subset_file <- file.path(data_dir, "pedigree_subset.csv")
dxfamily_file <- file.path(data_dir, "dxfamily.csv")

ped <- read.csv(ped_file, header = TRUE, stringsAsFactors = FALSE)
subset_ids_df <- read.csv(subset_file, header = TRUE, stringsAsFactors = FALSE)
dxfamily <- read.csv(dxfamily_file, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
if (!all(c("animal", "dam", "sire") %in% names(ped))) {
  stop("ped.csv must contain animal, dam, and sire columns.")
}
if (!"animal" %in% names(dxfamily)) {
  stop("dxfamily.csv must contain an animal column.")
}
if (!"animal" %in% names(subset_ids_df)) {
  stop("pedigree_subset.csv must contain an animal column.")
}

# Part 1: Clean the pedigree and include parents without individual rows.
ped_clean <- ped %>%
  transmute(
    animal = trimws(as.character(animal)),
    dam = trimws(as.character(dam)),
    sire = trimws(as.character(sire))
  ) %>%
  mutate(
    animal = na_if(animal, ""),
    dam = na_if(dam, ""),
    sire = na_if(sire, ""),
    dam = na_if(dam, "0"),
    sire = na_if(sire, "0"),
    dam = na_if(dam, "NA"),
    sire = na_if(sire, "NA")
  ) %>%
  filter(!is.na(animal), animal != "NA") %>%
  distinct(animal, .keep_all = TRUE)
all_ids <- unique(c(ped_clean$animal, ped_clean$dam, ped_clean$sire))
all_ids <- all_ids[!is.na(all_ids) & all_ids != "" & all_ids != "NA"]

ped_full <- tibble(animal = all_ids) %>%
  left_join(ped_clean, by = "animal") %>%
  distinct(animal, .keep_all = TRUE)

# Part 2: Calculate pedigree depth and equivalent complete generations.
ped_depth <- function(ped_df) {
  dam <- setNames(ped_df$dam, ped_df$animal)
  sire <- setNames(ped_df$sire, ped_df$animal)
  memo <- new.env(parent = emptyenv())
  
  f_depth <- function(id) {
    if (is.na(id) || !nzchar(id)) return(0L)
    if (exists(id, envir = memo, inherits = FALSE)) return(get(id, envir = memo))
    
    d <- dam[[id]]
    s <- sire[[id]]
    val <- 0L
    
    if (!is.null(d) && !is.na(d)) val <- max(val, 1L + f_depth(d))
    if (!is.null(s) && !is.na(s)) val <- max(val, 1L + f_depth(s))
    
    assign(id, val, envir = memo)
    val
  }
  
  tibble(
    animal = ped_df$animal,
    depth = vapply(ped_df$animal, f_depth, integer(1))
  )
}

ped_ecg <- function(ped_df) {
  dam <- setNames(ped_df$dam, ped_df$animal)
  sire <- setNames(ped_df$sire, ped_df$animal)
  memo <- new.env(parent = emptyenv())
  
  f_ecg <- function(id) {
    if (is.na(id) || !nzchar(id)) return(0)
    if (exists(id, envir = memo, inherits = FALSE)) return(get(id, envir = memo))
    
    d <- dam[[id]]
    s <- sire[[id]]
    val <- 0
    
    if (!is.null(d) && !is.na(d)) val <- val + (1 + f_ecg(d)) / 2
    if (!is.null(s) && !is.na(s)) val <- val + (1 + f_ecg(s)) / 2
    
    assign(id, val, envir = memo)
    val
  }
  
  tibble(
    animal = ped_df$animal,
    ecg = vapply(ped_df$animal, f_ecg, numeric(1))
  )
}

df_depth <- ped_depth(ped_full)
df_ecg <- ped_ecg(ped_full)

ped_metrics <- ped_full %>%
  select(animal, dam, sire) %>%
  left_join(df_depth, by = "animal") %>%
  left_join(df_ecg, by = "animal") %>%
  mutate(
    has_dam = !is.na(dam),
    has_sire = !is.na(sire),
    has_both = has_dam & has_sire,
    founder = !has_dam & !has_sire
  )

# Part 3: Calculate inbreeding coefficients from the full pedigree.
pedA <- ped_metrics %>%
  select(animal, dam, sire) %>%
  rename(id = animal) %>%
  as.data.frame()

pedA$dam[pedA$dam %in% c("0", "", "NA")] <- NA
pedA$sire[pedA$sire %in% c("0", "", "NA")] <- NA

invA <- MCMCglmm::inverseA(pedA, nodes = "ALL")
F_all <- setNames(invA$inbreeding, pedA$id)

# Part 4: Identify the pedigree subset and breeding females 
ids_pheno <- trimws(as.character(subset_ids_df$animal))
ids_pheno <- unique(ids_pheno[!is.na(ids_pheno) & ids_pheno != "" & ids_pheno != "NA"])

ids_female <- dxfamily %>%
  mutate(animal = trimws(as.character(animal))) %>%
  filter(!is.na(animal), animal != "", animal != "NA") %>%
  filter(if ("ld41" %in% names(.)) !is.na(ld41) else TRUE) %>%
  distinct(animal) %>%
  pull(animal)

ped_ph <- ped_metrics %>% filter(animal %in% ids_pheno)
ped_female <- ped_metrics %>% filter(animal %in% ids_female)

F_pheno <- F_all[names(F_all) %in% ids_pheno]
F_female <- F_all[names(F_all) %in% ids_female]
missing_pheno_in_ped <- setdiff(ids_pheno, ped_metrics$animal)
missing_female_in_ped <- setdiff(ids_female, ped_metrics$animal)

cat("Phenotyped pedigree subset IDs:", length(ids_pheno), "\n")
cat("Phenotyped pedigree subset IDs missing in pedigree:", length(missing_pheno_in_ped), "\n")
cat("Phenotyped females:", length(ids_female), "\n")
cat("Phenotyped females missing in pedigree:", length(missing_female_in_ped), "\n")

# Part 5: Summarize pedigree statistics.
summarise_pedigree_column <- function(ped_subset, F_subset) {
  c(
    nrow(ped_subset),
    sum(ped_subset$founder),
    sum(ped_subset$has_dam & !ped_subset$has_sire) +
      sum(!ped_subset$has_dam & ped_subset$has_sire),
    sum(ped_subset$has_both),
    round(mean(ped_subset$has_dam), 4),
    round(mean(ped_subset$has_sire), 4),
    round(mean(ped_subset$has_both), 4),
    max(ped_subset$depth),
    sprintf("%.2f ± %.2f", mean(ped_subset$depth), sd(ped_subset$depth)),
    sprintf("%.2f ± %.2f", mean(ped_subset$ecg), sd(ped_subset$ecg)),
    round(max(ped_subset$ecg), 2),
    sprintf("%.4f ± %.4f", mean(F_subset, na.rm = TRUE), sd(F_subset, na.rm = TRUE)),
    sprintf(
      "%.4f [%.4f, %.4f]",
      median(F_subset, na.rm = TRUE),
      quantile(F_subset, 0.25, na.rm = TRUE),
      quantile(F_subset, 0.75, na.rm = TRUE)
    ),
    round(max(F_subset, na.rm = TRUE), 4)
  )
}

table_s1_all_rows <- tibble(
  Parameter = c(
    "Total number of individuals",
    "Number of founders",
    "Number of individuals with one known parent",
    "Number of individuals with both parents known",
    "Proportion with known dam",
    "Proportion with known sire",
    "Proportion with both parents known",
    "Pedigree depth: maximum",
    "Pedigree depth: mean ± SD",
    "Equivalent complete generations: mean ± SD",
    "Equivalent complete generations: maximum",
    "Inbreeding coefficient: mean ± SD",
    "Inbreeding coefficient: median [IQR]",
    "Inbreeding coefficient: maximum"
  ),
  `Full pedigree` = summarise_pedigree_column(ped_metrics, F_all),
  `Phenotyped pedigree subset` = summarise_pedigree_column(ped_ph, F_pheno),
  `Phenotyped breeding females` = summarise_pedigree_column(ped_female, F_female)
)

print(table_s1_all_rows)

# Part 6: Select Table S1 rows.
table_s1_final <- table_s1_all_rows %>%
  filter(Parameter %in% c(
    "Total number of individuals",
    "Number of founders",
    "Number of individuals with both parents known",
    "Proportion with both parents known",
    "Pedigree depth: maximum",
    "Pedigree depth: mean ± SD",
    "Equivalent complete generations: mean ± SD",
    "Equivalent complete generations: maximum",
    "Inbreeding coefficient: mean ± SD",
    "Inbreeding coefficient: maximum"
  )) %>%
  mutate(
    Parameter = recode(
      Parameter,
      "Total number of individuals" = "No. of individuals",
      "Number of founders" = "No. of founders",
      "Number of individuals with both parents known" = "No. with both parents known",
      "Pedigree depth: maximum" = "Max. pedigree depth",
      "Pedigree depth: mean ± SD" = "Pedigree depth (mean ± SD)",
      "Equivalent complete generations: mean ± SD" = "Equivalent complete generations (mean ± SD)",
      "Equivalent complete generations: maximum" = "Max. equivalent complete generations",
      "Inbreeding coefficient: mean ± SD" = "Inbreeding coefficient (mean ± SD)",
      "Inbreeding coefficient: maximum" = "Max. inbreeding coefficient"
    )
  )

print(table_s1_final)
# Part 7: Export pedigree summaries.
write.csv(
  table_s1_all_rows,
  file.path(output_dir, "part5_tableS1_pedigree_summary_all_rows.csv"),
  row.names = FALSE
)

write.csv(
  table_s1_final,
  file.path(output_dir, "part5_tableS1_pedigree_summary_final.csv"),
  row.names = FALSE
)

# Part 8: Extract WorldClim monthly data for 1960-2019.
library(terra)

worldclim_dir <- file.path(public_dir, "external_data", "worldclim", "historical_monthly")

site <- terra::vect(
  data.frame(
    lon = 91.09861667,
    lat = 30.49945
  ),
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

extract_monthly_values <- function(variable) {
  files <- list.files(
    worldclim_dir,
    pattern = paste0(
      "_", variable,
      "_[0-9]{4}-[0-9]{2}\\.tif$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  file_names <- basename(files)
  
  year <- as.integer(
    sub(
      ".*_([0-9]{4})-([0-9]{2})\\.tif$",
      "\\1",
      file_names
    )
  )
  
  month <- as.integer(
    sub(
      ".*_([0-9]{4})-([0-9]{2})\\.tif$",
      "\\2",
      file_names
    )
  )
  file_info <- data.frame(
    file = files,
    year = year,
    month = month,
    stringsAsFactors = FALSE
  )
  file_info <- file_info[
    file_info$year >= 1960 &
      file_info$year <= 2019,
  ]
  file_info <- file_info[
    order(file_info$year, file_info$month),
  ]
  stopifnot(nrow(file_info) == 720)
  stopifnot(all(table(file_info$year) == 12))
  stopifnot(
    !anyDuplicated(
      paste(file_info$year, file_info$month)
    )
  )
  climate_rasters <- terra::rast(file_info$file)
  site_values <- terra::extract(
    climate_rasters,
    site
  )
  file_info$value <- as.numeric(site_values[1, -1])
  stopifnot(!anyNA(file_info$value))
  
  return(file_info[, c("year", "month", "value")])
}

tmin_monthly <- extract_monthly_values("tmin")
tmax_monthly <- extract_monthly_values("tmax")
prec_monthly <- extract_monthly_values("prec")

stopifnot(
  identical(
    tmin_monthly[, c("year", "month")],
    tmax_monthly[, c("year", "month")]
  )
)

stopifnot(
  identical(
    tmin_monthly[, c("year", "month")],
    prec_monthly[, c("year", "month")]
  )
)

monthly_climate <- data.frame(
  year = tmin_monthly$year,
  month = tmin_monthly$month,
  tmin = tmin_monthly$value,
  tmax = tmax_monthly$value,
  prec = prec_monthly$value
)

# Part 9: Calculate annual mean temperature and total precipitation.
monthly_climate$tmean <- (
  monthly_climate$tmin +
    monthly_climate$tmax
) / 2

annual_temperature <- aggregate(
  tmean ~ year,
  data = monthly_climate,
  FUN = mean
)

names(annual_temperature)[2] <- "annual_mean_temperature"

annual_precipitation <- aggregate(
  prec ~ year,
  data = monthly_climate,
  FUN = sum
)

names(annual_precipitation)[2] <- "annual_total_precipitation"

# Part 10: Test annual climate trends and estimate slopes with 95% confidence intervals.
temperature_cor_test <- cor.test(
  annual_temperature$year,
  annual_temperature$annual_mean_temperature,
  method = "pearson"
)

temperature_r <- as.numeric(
  temperature_cor_test$estimate
)

temperature_p <- temperature_cor_test$p.value

precipitation_cor_test <- cor.test(
  annual_precipitation$year,
  annual_precipitation$annual_total_precipitation,
  method = "pearson"
)

precipitation_r <- as.numeric(
  precipitation_cor_test$estimate
)

precipitation_p <- precipitation_cor_test$p.value

temperature_lm <- lm(
  annual_mean_temperature ~ year,
  data = annual_temperature
)

precipitation_lm <- lm(
  annual_total_precipitation ~ year,
  data = annual_precipitation
)

temperature_slope <- coef(temperature_lm)["year"]
temperature_slope_ci <- confint(temperature_lm)["year", ]

precipitation_slope <- coef(precipitation_lm)["year"]
precipitation_slope_ci <- confint(precipitation_lm)["year", ]

cat(
  "Annual mean temperature trend, 1960-2019:\n",
  "n =", nrow(annual_temperature), "years\n",
  "r =", round(temperature_r, 4), "\n",
  "P =", format.pval(temperature_p, digits = 4), "\n",
  "slope =", round(temperature_slope, 4), "degrees C/year\n",
  "95% CI =", round(temperature_slope_ci[1], 4),
  "to", round(temperature_slope_ci[2], 4), "\n\n"
)

cat(
  "Annual total precipitation trend, 1960-2019:\n",
  "n =", nrow(annual_precipitation), "years\n",
  "r =", round(precipitation_r, 4), "\n",
  "P =", round(precipitation_p, 4), "\n",
  "slope =", round(precipitation_slope, 4), "mm/year\n",
  "95% CI =", round(precipitation_slope_ci[1], 4),
  "to", round(precipitation_slope_ci[2], 4), "\n"
)
