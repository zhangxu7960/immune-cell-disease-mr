library(readr)
library(dplyr)
library(data.table)


project_path <- "E:/immune_disease"
immune_disease_path <- file.path(project_path, "mr_result/1e5ivw")
immune_disease_sensitivity_path <- file.path(project_path, "mr_result/1e5")
immune_inflammatory_path <- file.path(project_path, "91_inflammatory_factors/immune_to_inflammatory/1e5ivw")
immune_inflammatory_sensitivity_path <- file.path(project_path, "91_inflammatory_factors/immune_to_inflammatory/1e5")
inflammatory_disease_path <- file.path(project_path, "91_inflammatory_factors/inflammatory_to_disease/1e5ivw_all_disease")
inflammatory_disease_sensitivity_path <- file.path(project_path, "91_inflammatory_factors/inflammatory_to_disease/1e5_all_disease")
output_path <- file.path(project_path, "91_inflammatory_factors/result_filter")
input_paths <- c(immune_disease_path, immune_disease_sensitivity_path,
                 immune_inflammatory_path, immune_inflammatory_sensitivity_path,
                 inflammatory_disease_path, inflammatory_disease_sensitivity_path)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
required_columns <- c("id.exposure", "id.outcome", "exposure", "outcome", "method", "b", "se", "pval")
id_columns <- cols(id.exposure = col_character(), id.outcome = col_character())
empty_result <- data.frame(id.exposure = character(), id.outcome = character(),
                           exposure = character(), outcome = character(), method = character(),
                           b = numeric(), se = numeric(), pval = numeric(),
                           pleiotropy_pval = numeric(), Q_pval = numeric())

# Immune cell -> disease
ivw_path <- immune_disease_path
sensitivity_path <- immune_disease_sensitivity_path
ivw_file <- list.files(ivw_path, pattern = "\\.csv$", full.names = TRUE)
ivw_result <- list(empty_result)
seen_tests <- list()
for (i in ivw_file) {
  data <- read_csv(i, col_types = id_columns, show_col_types = FALSE)
  data <- data[data$method %in% "Inverse variance weighted", ]
  if(nrow(data) == 0) next
  if(anyNA(data$id.exposure) || anyNA(data$id.outcome) ||
     any(!nzchar(trimws(data$id.exposure))) || any(!nzchar(trimws(data$id.outcome)))) stop("Missing GWAS IDs in: ", i)
  seen_tests[[length(seen_tests) + 1]] <- data %>% dplyr::select(id.exposure, id.outcome)
  data_clean <- subset(data, is.finite(b) & is.finite(pval) & pval >= 0 & pval < 0.01)
  if(nrow(data_clean) == 0) next
  outcome_id <- data_clean$id.outcome[1]
  exposure_ids <- unique(data_clean$id.exposure)
  all_ple <- list()
  all_heterogeneity <- list()
  for (exposure_id in exposure_ids) {
    ple_file <- file.path(sensitivity_path, outcome_id, exposure_id, "ple.csv")
    heterogeneity_file <- file.path(sensitivity_path, outcome_id, exposure_id, "H.csv")
    if(!file.exists(ple_file) || !file.exists(heterogeneity_file)) {
      message("Sensitivity results unavailable; excluded: ", outcome_id, "/", exposure_id)
      next
    }
    ple <- read_csv(ple_file, col_types = id_columns, show_col_types = FALSE)
    heterogeneity <- read_csv(heterogeneity_file, col_types = id_columns, show_col_types = FALSE)
    ple <- ple %>% dplyr::filter(id.exposure == exposure_id, id.outcome == outcome_id)
    heterogeneity <- heterogeneity %>% dplyr::filter(id.exposure == exposure_id, id.outcome == outcome_id, method == "MR Egger")
    if(nrow(ple) == 0 || nrow(heterogeneity) == 0) {
      message("Required sensitivity tests unavailable; excluded: ", outcome_id, "/", exposure_id)
      next
    }
    all_ple[[exposure_id]] <- ple %>%
      dplyr::select(id.exposure, id.outcome, pval) %>% dplyr::rename(pleiotropy_pval = pval)
    all_heterogeneity[[exposure_id]] <- heterogeneity %>%
      dplyr::select(id.exposure, id.outcome, Q_pval)
  }
  if(length(all_ple) == 0) next
  ple_combined <- dplyr::bind_rows(all_ple)
  heterogeneity_combined <- dplyr::bind_rows(all_heterogeneity)
  ple_heterogeneity_combined <- dplyr::inner_join(ple_combined, heterogeneity_combined, by = c("id.exposure", "id.outcome"))
  data_clean <- data_clean %>% dplyr::select(-dplyr::any_of(c("pleiotropy_pval", "Q_pval")))
  result <- dplyr::inner_join(data_clean, ple_heterogeneity_combined, by = c("id.exposure", "id.outcome"))
  result_clean <- subset(result, is.finite(pleiotropy_pval) & pleiotropy_pval > 0.05 & pleiotropy_pval <= 1 &
                                  is.finite(Q_pval) & Q_pval > 0.05 & Q_pval <= 1)
  if(nrow(result_clean) > 0) ivw_result[[length(ivw_result) + 1]] <- result_clean
}
immune_disease_result <- dplyr::bind_rows(ivw_result)
fwrite(immune_disease_result, file = file.path(output_path, "immune_disease_p_0.01.csv"), row.names = FALSE)
message("Immune cell -> disease: ", nrow(immune_disease_result), " associations retained.")

# Immune cell -> inflammatory
ivw_path <- immune_inflammatory_path
sensitivity_path <- immune_inflammatory_sensitivity_path
ivw_file <- list.files(ivw_path, pattern = "\\.csv$", full.names = TRUE)
ivw_result <- list(empty_result)
seen_tests <- list()
for (i in ivw_file) {
  data <- read_csv(i, col_types = id_columns, show_col_types = FALSE)
  if(!all(required_columns %in% names(data))) stop("Missing MR columns in: ", i)
  data <- data[data$method %in% "Inverse variance weighted", ]
  if(nrow(data) == 0) next
  seen_tests[[length(seen_tests) + 1]] <- data %>% dplyr::select(id.exposure, id.outcome)
  data_clean <- subset(data, is.finite(b) & is.finite(pval) & pval >= 0 & pval < 0.01)
  if(nrow(data_clean) == 0) next
  outcome_id <- data_clean$id.outcome[1]
  exposure_ids <- unique(data_clean$id.exposure)
  all_ple <- list()
  all_heterogeneity <- list()
  for (exposure_id in exposure_ids) {
    ple_file <- file.path(sensitivity_path, outcome_id, exposure_id, "ple.csv")
    heterogeneity_file <- file.path(sensitivity_path, outcome_id, exposure_id, "H.csv")
    if(!file.exists(ple_file) || !file.exists(heterogeneity_file)) {
      message("Sensitivity results unavailable; excluded: ", outcome_id, "/", exposure_id)
      next
    }
    ple <- read_csv(ple_file, col_types = id_columns, show_col_types = FALSE)
    heterogeneity <- read_csv(heterogeneity_file, col_types = id_columns, show_col_types = FALSE)
    ple <- ple %>% dplyr::filter(id.exposure == exposure_id, id.outcome == outcome_id)
    heterogeneity <- heterogeneity %>% dplyr::filter(id.exposure == exposure_id, id.outcome == outcome_id, method == "MR Egger")
    if(nrow(ple) == 0 || nrow(heterogeneity) == 0) {
      message("Required sensitivity tests unavailable; excluded: ", outcome_id, "/", exposure_id)
      next
    }
    all_ple[[exposure_id]] <- ple %>%
      dplyr::select(id.exposure, id.outcome, pval) %>% dplyr::rename(pleiotropy_pval = pval)
    all_heterogeneity[[exposure_id]] <- heterogeneity %>%
      dplyr::select(id.exposure, id.outcome, Q_pval)
  }
  if(length(all_ple) == 0) next
  ple_combined <- dplyr::bind_rows(all_ple)
  heterogeneity_combined <- dplyr::bind_rows(all_heterogeneity)
  ple_heterogeneity_combined <- dplyr::inner_join(ple_combined, heterogeneity_combined, by = c("id.exposure", "id.outcome"))
  data_clean <- data_clean %>% dplyr::select(-dplyr::any_of(c("pleiotropy_pval", "Q_pval")))
  result <- dplyr::inner_join(data_clean, ple_heterogeneity_combined, by = c("id.exposure", "id.outcome"))
  result_clean <- subset(result, is.finite(pleiotropy_pval) & pleiotropy_pval > 0.05 & pleiotropy_pval <= 1 &
                                  is.finite(Q_pval) & Q_pval > 0.05 & Q_pval <= 1)
  if(nrow(result_clean) > 0) ivw_result[[length(ivw_result) + 1]] <- result_clean
}
immune_inflammatory_result <- dplyr::bind_rows(ivw_result)
fwrite(immune_inflammatory_result, file = file.path(output_path, "immune_inflammatory_p_0.01.csv"), row.names = FALSE)
message("Immune cell -> inflammatory: ", nrow(immune_inflammatory_result), " associations retained.")

# Inflammatory -> disease
ivw_path <- inflammatory_disease_path
sensitivity_path <- inflammatory_disease_sensitivity_path
ivw_file <- list.files(ivw_path, pattern = "\\.csv$", full.names = TRUE)
if(length(ivw_file) == 0) stop("No MR summary CSV files found: ", ivw_path)
ivw_result <- list(empty_result)
seen_tests <- list()
for (i in ivw_file) {
  data <- read_csv(i, col_types = id_columns, show_col_types = FALSE)
  data <- data[data$method %in% "Inverse variance weighted", ]
  if(nrow(data) == 0) next
  seen_tests[[length(seen_tests) + 1]] <- data %>% dplyr::select(id.exposure, id.outcome)
  data_clean <- subset(data, is.finite(b) & is.finite(pval) & pval >= 0 & pval < 0.01)
  if(nrow(data_clean) == 0) next
  outcome_id <- data_clean$id.outcome[1]
  exposure_ids <- unique(data_clean$id.exposure)
  all_ple <- list()
  all_heterogeneity <- list()
  for (exposure_id in exposure_ids) {
    ple_file <- file.path(sensitivity_path, outcome_id, exposure_id, "ple.csv")
    heterogeneity_file <- file.path(sensitivity_path, outcome_id, exposure_id, "H.csv")
    if(!file.exists(ple_file) || !file.exists(heterogeneity_file)) {
      message("Sensitivity results unavailable; excluded: ", outcome_id, "/", exposure_id)
      next
    }
    ple <- read_csv(ple_file, col_types = id_columns, show_col_types = FALSE)
    heterogeneity <- read_csv(heterogeneity_file, col_types = id_columns, show_col_types = FALSE)
    ple <- ple %>% dplyr::filter(id.exposure == exposure_id, id.outcome == outcome_id)
    heterogeneity <- heterogeneity %>% dplyr::filter(id.exposure == exposure_id, id.outcome == outcome_id, method == "MR Egger")
    if(nrow(ple) == 0 || nrow(heterogeneity) == 0) {
      message("Required sensitivity tests unavailable; excluded: ", outcome_id, "/", exposure_id)
      next
    }
    all_ple[[exposure_id]] <- ple %>%
      dplyr::select(id.exposure, id.outcome, pval) %>% dplyr::rename(pleiotropy_pval = pval)
    all_heterogeneity[[exposure_id]] <- heterogeneity %>%
      dplyr::select(id.exposure, id.outcome, Q_pval)
  }
  if(length(all_ple) == 0) next
  ple_combined <- dplyr::bind_rows(all_ple)
  heterogeneity_combined <- dplyr::bind_rows(all_heterogeneity)
  ple_heterogeneity_combined <- dplyr::inner_join(ple_combined, heterogeneity_combined, by = c("id.exposure", "id.outcome"))
  data_clean <- data_clean %>% dplyr::select(-dplyr::any_of(c("pleiotropy_pval", "Q_pval")))
  result <- dplyr::inner_join(data_clean, ple_heterogeneity_combined, by = c("id.exposure", "id.outcome"))
  result_clean <- subset(result, is.finite(pleiotropy_pval) & pleiotropy_pval > 0.05 & pleiotropy_pval <= 1 &
                                  is.finite(Q_pval) & Q_pval > 0.05 & Q_pval <= 1)
  if(nrow(result_clean) > 0) ivw_result[[length(ivw_result) + 1]] <- result_clean
}
inflammatory_disease_result <- dplyr::bind_rows(ivw_result)
fwrite(inflammatory_disease_result, file = file.path(output_path, "inflammatory_disease_p_0.01.csv"), row.names = FALSE)
message("Inflammatory -> disease: ", nrow(inflammatory_disease_result), " associations retained.")

# Match the three associations by GWAS IDs and retain the sensitivity P values for each association.
immune_inflammatory <- immune_inflammatory_result %>%
  dplyr::select(id.exposure, id.outcome, outcome, b, se, pval, pleiotropy_pval, Q_pval) %>%
  dplyr::rename(id.immune_cell = id.exposure, id.inflammatory = id.outcome, inflammatory = outcome,
                b1 = b, se1 = se, pval1 = pval, pleiotropy_pval1 = pleiotropy_pval, Q_pval1 = Q_pval)
inflammatory_disease <- inflammatory_disease_result %>%
  dplyr::select(id.exposure, id.outcome, b, se, pval, pleiotropy_pval, Q_pval) %>%
  dplyr::rename(id.inflammatory = id.exposure, id.disease = id.outcome,
                b2 = b, se2 = se, pval2 = pval, pleiotropy_pval2 = pleiotropy_pval, Q_pval2 = Q_pval)
immune_disease <- immune_disease_result %>%
  dplyr::select(id.exposure, exposure, id.outcome, outcome, b, se, pval, pleiotropy_pval, Q_pval,
                dplyr::any_of(c("cell_type", "category"))) %>%
  dplyr::rename(id.immune_cell = id.exposure, immune_cell = exposure, id.disease = id.outcome, disease = outcome)
immune_inflammatory_disease <- inner_join(immune_inflammatory, inflammatory_disease, by = "id.inflammatory")
media_result <- inner_join(immune_inflammatory_disease, immune_disease, by = c("id.immune_cell", "id.disease")) %>%
  mutate(mediated_effect = b1 * b2)

# Retain same-direction indirect and total effects and a mediation fraction strictly between 0 and 1.
media_mr_final_result <- media_result %>%
  filter(is.finite(mediated_effect), (mediated_effect > 0 & b > 0) | (mediated_effect < 0 & b < 0)) %>%
  mutate(direct_effect = b - mediated_effect, mediated_proportion = mediated_effect / b) %>%
  filter(is.finite(mediated_proportion), mediated_proportion > 0, mediated_proportion < 1) %>%
  dplyr::select(id.immune_cell, immune_cell, id.inflammatory, inflammatory, id.disease, disease,
                b1, se1, pval1, b2, se2, pval2, b, se, pval,
                mediated_effect, direct_effect, mediated_proportion,
                pleiotropy_pval1, Q_pval1, pleiotropy_pval2, Q_pval2, pleiotropy_pval, Q_pval,
                dplyr::any_of(c("cell_type", "category"))) %>%
  arrange(id.immune_cell, id.inflammatory, id.disease)
fwrite(media_mr_final_result, file = file.path(output_path, "inflammatory_final_result.csv"), row.names = FALSE)

