library(dplyr)
library(readr)


data_path <- "E:/immune_disease/sceQTL_MR_14_immune_cell_types/1e5ivw_0.001_10000kb_clean"
result_path <- "E:/immune_disease/mental_disorders_B_cells/eGene_cross_disease_associations"
outcome_id <- c("DGWAS-166", "DGWAS-2628", "DGWAS-6186",
                "DGWAS-481", "DGWAS-4147")
plot_columns <- c("id.outcome", "exposure", "cell_type", "b", "FDR")
cell_order <- c("bin", "bmem", "cd4et", "cd4nc", "cd4sox4", "cd8et", "cd8nc",
                "cd8s100b", "dc", "monoc", "mononc", "nk", "nkr", "plasma")


mr_list <- lapply(outcome_id, function(id) {
  dat <- read_csv(file.path(data_path, paste0(id, ".csv")), show_col_types = FALSE)
  dat %>% dplyr::select(all_of(plot_columns))
})
plot_result <- bind_rows(mr_list)
plot_result <- plot_result %>%
  arrange(match(id.outcome, outcome_id), match(cell_type, cell_order), exposure)


dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
write_excel_csv(plot_result, file.path(result_path, "eGene_plot_data.csv"))


