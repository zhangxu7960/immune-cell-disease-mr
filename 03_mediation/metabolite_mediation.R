library(readr)
library(dplyr)
library(data.table)

project_path <- "E:/immune_disease"
immune_disease_path <- file.path(project_path, "mr_result/1e5ivw")
immune_disease_sensitivity_path <- file.path(project_path, "mr_result/1e5")
immune_metabolite_path <- file.path(project_path, "1400_daixiewu/immune_metabolite/mr_result/1e5ivw_raps")
immune_metabolite_sensitivity_path <- file.path(project_path, "1400_daixiewu/immune_metabolite/mr_result/1e5_raps")
metabolite_disease_path <- file.path(project_path, "1400_daixiewu/metabolite_disease/mr_result/1e5ivw")
metabolite_disease_sensitivity_path <- file.path(project_path, "1400_daixiewu/metabolite_disease/mr_result/1e5")
output_path <- file.path(project_path, "1400_daixiewu/filter_result_p_0.01/result_1")
input_paths <- c(immune_disease_path, immune_disease_sensitivity_path,
                 immune_metabolite_path, immune_metabolite_sensitivity_path,
                 metabolite_disease_path, metabolite_disease_sensitivity_path)
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

# Immune cell -> metabolite
ivw_path <- immune_metabolite_path
sensitivity_path <- immune_metabolite_sensitivity_path
ivw_file <- list.files(ivw_path, pattern = "\\.csv$", full.names = TRUE)
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
immune_metabolite_result <- dplyr::bind_rows(ivw_result)
fwrite(immune_metabolite_result, file = file.path(output_path, "immune_metabolite_p_0.01.csv"), row.names = FALSE)
message("Immune cell -> metabolite: ", nrow(immune_metabolite_result), " associations retained.")

# Metabolite -> disease
ivw_path <- metabolite_disease_path
sensitivity_path <- metabolite_disease_sensitivity_path
ivw_file <- list.files(ivw_path, pattern = "\\.csv$", full.names = TRUE)
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
    if(nrow(ple) > 1 || nrow(heterogeneity) > 1) stop("Duplicate sensitivity results for: ", outcome_id, "/", exposure_id)
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
metabolite_disease_result <- dplyr::bind_rows(ivw_result)
fwrite(metabolite_disease_result, file = file.path(output_path, "metabolite_disease_p_0.01.csv"), row.names = FALSE)
message("Metabolite -> disease: ", nrow(metabolite_disease_result), " associations retained.")

# Match the three associations by GWAS IDs and retain the sensitivity P values for each association.
immune_metabolite <- immune_metabolite_result %>%
  dplyr::select(id.exposure, id.outcome, outcome, b, se, pval, pleiotropy_pval, Q_pval) %>%
  dplyr::rename(id.immune_cell = id.exposure, id.metabolite = id.outcome, metabolite = outcome,
                b1 = b, se1 = se, pval1 = pval, pleiotropy_pval1 = pleiotropy_pval, Q_pval1 = Q_pval)
metabolite_disease <- metabolite_disease_result %>%
  dplyr::select(id.exposure, id.outcome, b, se, pval, pleiotropy_pval, Q_pval) %>%
  dplyr::rename(id.metabolite = id.exposure, id.disease = id.outcome,
                b2 = b, se2 = se, pval2 = pval, pleiotropy_pval2 = pleiotropy_pval, Q_pval2 = Q_pval)
immune_disease <- immune_disease_result %>%
  dplyr::select(id.exposure, exposure, id.outcome, outcome, b, se, pval, pleiotropy_pval, Q_pval,
                dplyr::any_of(c("cell_type", "category"))) %>%
  dplyr::rename(id.immune_cell = id.exposure, immune_cell = exposure, id.disease = id.outcome, disease = outcome)
immune_metabolite_disease <- inner_join(immune_metabolite, metabolite_disease, by = "id.metabolite")
media_result <- inner_join(immune_metabolite_disease, immune_disease, by = c("id.immune_cell", "id.disease")) %>%
  mutate(mediated_effect = b1 * b2)

media_mr_final_result <- media_result %>%
  filter(is.finite(mediated_effect), (mediated_effect > 0 & b > 0) | (mediated_effect < 0 & b < 0)) %>%
  mutate(direct_effect = b - mediated_effect, mediated_proportion = mediated_effect / b) %>%
  filter(is.finite(mediated_proportion), mediated_proportion > 0, mediated_proportion < 1) %>%
  dplyr::select(id.immune_cell, immune_cell, id.metabolite, metabolite, id.disease, disease,
                b1, se1, pval1, b2, se2, pval2, b, se, pval,
                mediated_effect, direct_effect, mediated_proportion,
                pleiotropy_pval1, Q_pval1, pleiotropy_pval2, Q_pval2, pleiotropy_pval, Q_pval,
                dplyr::any_of(c("cell_type", "category"))) %>%
  arrange(id.immune_cell, id.metabolite, id.disease)
fwrite(media_mr_final_result, file = file.path(output_path, "media_result_final.csv"), row.names = FALSE)

