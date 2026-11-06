library(dplyr)
library(readr)


mr_file <- "E:/immune_disease/dementia_monocytes/all_cells_significant_eGene_associations.csv"
result_path <- "E:/immune_disease/dementia_monocytes/candidate_gene_MR_forest"
result_file <- file.path(result_path, "Monocyte_candidate_gene_MR_associations.csv")

target_cells <- c("monoc", "mononc")
cell_labels <- c(monoc = "MONO C", mononc = "MONO NC")
outcome_id <- "DGWAS-105"
outcome_labels <- "Dementia"
fdr_cutoff <- 0.05
ci_multiplier <- 1.96


mr_result <- read_csv(mr_file, show_col_types = FALSE)
required_columns <- c("exposure", "cell_type", "id.outcome", "method", "nsnp",
                      "b", "se", "or", "FDR")
mr_result <- mr_result %>% mutate(source_row = row_number() + 1L)
plot_result <- mr_result %>%
  dplyr::filter(cell_type %in% target_cells,
                id.outcome %in% outcome_id,
                is.finite(FDR), FDR < fdr_cutoff)


plot_result <- plot_result %>%
  mutate(outcome_order = match(id.outcome, outcome_id),
         disease_label_en = outcome_labels[outcome_order],
         cell_label_plot = unname(cell_labels[cell_type]),
         cell_order = match(cell_type, target_cells),
         or_lci95 = exp(b - ci_multiplier * se),
         or_uci95 = exp(b + ci_multiplier * se)) %>%
  arrange(outcome_order, exposure, cell_order) %>%
  mutate(plot_order = row_number()) %>%
  dplyr::select(-outcome_order, -cell_order)



dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
write_excel_csv(plot_result, result_file, na = "")

print(plot_result %>% count(disease_label_en, name = "n_associations"))


