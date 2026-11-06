library(readr)
library(dplyr)

# Within each outcome, apply BH separately to Wald ratio and IVW tests across all cell types.

base_path <- "E:/immune_disease"
ivw_path <- file.path(base_path, "1e5ivw_0.001_10000kb")
output_path <- file.path(base_path, "1e5ivw_0.001_10000kb_clean")
output_positive_path <- file.path(base_path, "1e5ivw_0.001_10000kb_positve")

files_path <- list.files(ivw_path, pattern = "\\.csv$", full.names = TRUE,
                         ignore.case = TRUE)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_positive_path, recursive = TRUE, showWarnings = FALSE)

result_all <- vector("list", length(files_path))
no_positive_result <- character()

for (i in seq_along(files_path)) {
  file_path <- files_path[[i]]
  data <- read_csv(file_path, show_col_types = FALSE)

  required_columns <- c("id.exposure", "id.outcome", "cell_type", "method", "pval")
  if (!all(required_columns %in% names(data))) {
    stop("Missing required columns in: ", file_path)
  }
  if (length(unique(data$id.outcome)) > 1L) {
    stop("Each input CSV must contain only one outcome: ", file_path)
  }

  data_clean <- data %>%
    dplyr::filter(method %in% c("Wald ratio", "Inverse variance weighted"),
                  is.finite(pval), pval >= 0, pval <= 1)
  if (anyDuplicated(data_clean[c("id.exposure", "id.outcome", "cell_type")]) > 0L) {
    stop("Duplicate primary gene-cell-outcome tests in: ", file_path)
  }
  data_clean <- data_clean %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(FDR = p.adjust(pval, method = "BH")) %>%
    dplyr::ungroup()
  write_csv(data_clean, file = file.path(output_path, basename(file_path)))

  positive_result <- data_clean %>% dplyr::filter(FDR < 0.05)
  # Write header-only files when no associations are significant to avoid stale results.
  write_csv(positive_result,
            file = file.path(output_positive_path, basename(file_path)))
  result_all[[i]] <- positive_result
  if (nrow(positive_result) == 0L) {
    no_positive_result <- c(no_positive_result, basename(file_path))
  }
  message(basename(file_path), ": ", nrow(data_clean), " tests; ",
          nrow(positive_result), " with FDR < 0.05")
}

result <- dplyr::bind_rows(result_all)
write_csv(result, file = file.path(base_path, "sceqtl_disease_fdr_0.05.csv"))
writeLines(no_positive_result,
           con = file.path(base_path, "no_positive_result.txt"), useBytes = TRUE)
message("Completed: ", length(files_path), " outcome files; ", nrow(result),
        " associations with FDR < 0.05.")
