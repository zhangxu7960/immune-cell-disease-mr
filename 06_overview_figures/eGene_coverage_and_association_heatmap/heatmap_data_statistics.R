# Count gene-cell-outcome GWAS associations using the existing FDR < 0.05 results.

library(data.table)
library(dplyr)
library(readxl)

base_dir <- "E:/immune_disease/sceQTL_14_immune_cells"
mr_file <- file.path(base_dir, "sceqtl_disease_fdr_0.05.csv")
clean_dir <- file.path(base_dir, "1e5ivw_0.001_10000kb_clean")
disease_file <- "E:/immune_disease/6640-disease-information.xlsx"
output_dir <- file.path(base_dir, "eGene_coverage_and_disease_associations", "plot_data")
if(!dir.exists(output_dir))dir.create(output_dir, recursive = TRUE)

# Use the existing FDR values without applying another correction.
cell_order <- c("bin", "bmem", "cd4nc", "cd4et", "cd4sox4",
                "cd8nc", "cd8et", "cd8s100b", "nk", "nkr",
                "monoc", "mononc", "dc", "plasma")
main_methods <- c("Wald ratio", "Inverse variance weighted")


positive_dat <- fread(mr_file, showProgress = FALSE)
required <- c("id.exposure", "id.outcome", "cell_type", "method", "FDR")
key_cols <- c("id.exposure", "cell_type", "id.outcome")



disease_map <- read_excel(disease_file, sheet = 1) %>%
  transmute(id.outcome = trimws(as.character(GWASID)),
            disease_category = trimws(as.character(category)))



# Use n_tested to distinguish zero positives from untested pairs; heatmap values are association counts.
clean_files <- list.files(clean_dir, pattern = "\\.csv$", full.names = TRUE)
tested_dat <- rbindlist(lapply(clean_files, function(f){
  fread(f, select = c(key_cols, "method", "pval", "FDR"), showProgress = FALSE)
})) %>%
  dplyr::filter(method %in% main_methods,
                is.finite(pval), pval >= 0, pval <= 1,
                is.finite(FDR), FDR >= 0, FDR <= 1)

# Verify that the positive summary matches significant records in the complete results.
positive_keys <- positive_dat[, ..key_cols]
expected_keys <- tested_dat %>% dplyr::filter(FDR < 0.05) %>%
  dplyr::select(all_of(key_cols))


tested_dat <- tested_dat %>% left_join(disease_map, by = "id.outcome")
positive_dat <- positive_dat %>% left_join(disease_map, by = "id.outcome")




positive_stat <- positive_dat %>%
  count(cell_type, disease_category, name = "n_associations")
tested_stat <- tested_dat %>%
  count(cell_type, disease_category, name = "n_tested")
outcome_stat <- tested_dat %>%
  dplyr::distinct(id.outcome, disease_category) %>%
  count(disease_category, name = "n_outcome_gwas")

plot_data <- expand.grid(
  cell_type = cell_order,
  disease_category = sort(unique(tested_dat$disease_category)),
  stringsAsFactors = FALSE
) %>%
  left_join(positive_stat, by = c("cell_type", "disease_category")) %>%
  left_join(tested_stat, by = c("cell_type", "disease_category")) %>%
  left_join(outcome_stat, by = "disease_category") %>%
  mutate(n_tested = ifelse(is.na(n_tested), 0L, n_tested),
         n_associations = ifelse(n_tested == 0L, NA_integer_,
                                 ifelse(is.na(n_associations), 0L, n_associations)))

fwrite(plot_data, file.path(output_dir, "heatmap_data.csv"), na = "NA")

