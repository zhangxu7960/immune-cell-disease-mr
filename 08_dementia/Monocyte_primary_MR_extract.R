library(dplyr)
library(readr)

mr_file <- "E:/immune_disease/immune_disease_MR/immune_disease_fdr_0.05_ple_het_reverse_20261012.csv"
result_path <- "E:/immune_disease/dementia_monocytes/primary_MR_forest"
result_file <- file.path(result_path, "Monocyte_dementia_primary_MR_associations.csv")

target_cell <- "Monocyte"
target_category <- "Mental and behavioural disorders"
outcome_id <- "DGWAS-105"
outcome_label <- "Dementia"
fdr_cutoff <- 0.05
ci_multiplier <- 1.96

phenotype_id <- c("GCST90001988", "GCST90001991", "GCST90002010")


mr_result <- read_csv(mr_file, show_col_types = FALSE)
required_columns <- c("id.exposure", "id.outcome", "exposure", "outcome",
                      "cell_type", "category", "method", "nsnp",
                      "b", "se", "or", "FDR")

mr_result <- mr_result %>% mutate(source_row = row_number() + 1L)



plot_result <- mr_result %>%
  dplyr::filter(id.exposure %in% phenotype_id, id.outcome == outcome_id)


plot_result <- plot_result %>%
  mutate(phenotype_order = match(id.exposure, phenotype_id),
         disease_label = outcome_label,
         or_lci95 = exp(b - ci_multiplier * se),
         or_uci95 = exp(b + ci_multiplier * se)) %>%
  arrange(phenotype_order) %>%
  mutate(plot_order = row_number()) %>%
  dplyr::select(-phenotype_order)



dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

write_excel_csv(plot_result, result_file, na = "")

print(plot_result %>% dplyr::select(exposure, disease_label, or,
                                   or_lci95, or_uci95, FDR))

