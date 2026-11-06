library(dplyr)
library(readr)


mr_file <- "E:/immune_disease/immune_disease_MR/immune_disease_fdr_0.05_ple_het_reverse_20261012.csv"
result_path <- "E:/immune_disease/mental_disorders_B_cells/primary_MR_forest"
sceqtl_path <- "E:/immune_disease/sceQTL_MR_14_immune_cell_types/1e5ivw_0.001_10000kb_clean"
main_result_file <- file.path(result_path,
  "B_cell_primary_MR_main_DGWAS-166_2628_6186_481.csv")
supp_result_file <- file.path(result_path,
  "B_cell_primary_MR_supplement_DGWAS-4147.csv")

target_cell <- "B cell"
target_category <- "Mental and behavioural disorders"
fdr_cutoff <- 0.05
ci_multiplier <- 1.96


outcome_id <- c("DGWAS-166", "DGWAS-2628", "DGWAS-6186",
                "DGWAS-481", "DGWAS-4147")
outcome_label <- c("Anxiety disorders", "Mental disorder",
                   "Depression-related phenotype", "Vascular dementia",
                   "Vascular dementia")


mr_result <- read_csv(mr_file, show_col_types = FALSE)
required_columns <- c("id.exposure", "id.outcome", "exposure", "outcome",
                      "cell_type", "category", "method", "nsnp",
                      "b", "se", "or", "FDR")

mr_result <- mr_result %>% mutate(source_row = row_number() + 1L)



plot_result <- mr_result %>%
  dplyr::filter(cell_type == target_cell,
                category == target_category,
                is.finite(FDR), FDR < fdr_cutoff)

plot_result <- plot_result %>%
  mutate(outcome_order = match(id.outcome, outcome_id),
         disease_label = outcome_label[outcome_order],
         or_lci95 = exp(b - ci_multiplier * se),
         or_uci95 = exp(b + ci_multiplier * se)) %>%
  arrange(outcome_order, FDR, id.exposure) %>%
  mutate(plot_order = row_number()) %>%
  dplyr::select(-outcome_order)





# Use the existing FDR values from the source sceQTL MR results.
# Count B IN and B MEM gene associations separately; count unique genes independently.
egene_evidence <- bind_rows(lapply(outcome_id, function(id) {
  sceqtl_file <- file.path(sceqtl_path, paste0(id, ".csv"))
  dat <- read_csv(sceqtl_file, show_col_types = FALSE)
  b_hits <- dat %>%
    filter(cell_type %in% c("bin", "bmem"),
           is.finite(FDR), FDR < fdr_cutoff)
  data.frame(id.outcome = id,
             B_cell_eGene_MR_association_n = nrow(b_hits),
             B_cell_eGene_n = n_distinct(b_hits$exposure))
}))

plot_result <- plot_result %>%
  left_join(egene_evidence, by = "id.outcome") %>%
  mutate(figure_role = if_else(B_cell_eGene_MR_association_n > 0L,
                               "Main figure", "Supplementary figure"))
main_result <- plot_result %>%
  filter(figure_role == "Main figure") %>% mutate(plot_order = row_number())
supp_result <- plot_result %>%
  filter(figure_role == "Supplementary figure") %>% mutate(plot_order = row_number())


stopifnot(nrow(main_result) == 9L, nrow(supp_result) == 1L,
          setequal(main_result$id.outcome,
                   c("DGWAS-166", "DGWAS-2628", "DGWAS-6186", "DGWAS-481")),
          identical(unique(supp_result$id.outcome), "DGWAS-4147"),
          nrow(main_result) + nrow(supp_result) == nrow(plot_result))

dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

write_excel_csv(main_result, main_result_file, na = "")
write_excel_csv(supp_result, supp_result_file, na = "")


