# Count unique exposure GWAS ID-outcome GWAS ID pairs.
# N: all eligible test pairs within a cell type.
# K: eligible test pairs in the disease system within that cell type.
# n: retained positive pairs within that cell type.
# k: retained positive pairs in the disease system within that cell type.
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(readxl)


immune_type <- read_xlsx("E:/immune_disease/731-immune-information.xlsx") %>%
  dplyr::rename(id.exposure = id) %>%
  dplyr::select(id.exposure, Cell_type)

disease_type <- read_xlsx("E:/immune_disease/6640-disease-information.xlsx") %>%
  dplyr::rename(id.outcome = GWASID) %>%
  dplyr::select(id.outcome, category) %>%
  filter(!is.na(category), category != "")
stopifnot(!anyDuplicated(immune_type$id.exposure), !anyDuplicated(disease_type$id.outcome))

# Apply the same classification exclusions to foreground and background pairs.
result <- read_csv("E:/immune_disease/immune_disease_fdr_0.05_ple_het_reverse.csv",
                   show_col_types = FALSE) %>%
  dplyr::select(id.exposure, id.outcome) %>%
  distinct(id.exposure, id.outcome) %>%
  inner_join(immune_type, by = "id.exposure") %>%
  inner_join(disease_type, by = "id.outcome") %>%
  filter(!is.na(Cell_type), Cell_type != "")


# Build the background from all eligible tests, including non-significant results.
mr_files <- list.files("E:/immune_disease/mr_result/1e5ivw_adjP",
                       pattern = "\\.csv$", full.names = TRUE)
mr_files <- mr_files[sub("\\.csv$", "", basename(mr_files)) %in% disease_type$id.outcome]
stopifnot(length(mr_files) > 0)
background_list <- result_list <- vector("list", length(mr_files))

for (i in seq_along(mr_files)) {
  background <- read_csv(mr_files[i], show_col_types = FALSE, progress = FALSE,
                         num_threads = 1, col_types = cols_only(
                           id.exposure = col_character(), id.outcome = col_character(),
                           method = col_character(), nsnp = col_double(),
                           b = col_double(), se = col_double(), pval = col_double())) %>%
    filter(method == "Inverse variance weighted",
           is.finite(nsnp), nsnp >= 2, is.finite(b), is.finite(se), se > 0,
           is.finite(pval), pval >= 0, pval <= 1) %>%
    distinct(id.exposure, id.outcome, .keep_all = TRUE)
  stopifnot(all(background$id.outcome == sub("\\.csv$", "", basename(mr_files[i]))))
  result_list[[i]] <- result %>% semi_join(background, by = c("id.exposure", "id.outcome"))
  background_list[[i]] <- background %>%
    inner_join(immune_type, by = "id.exposure") %>%
    filter(!is.na(Cell_type), Cell_type != "") %>%
    group_by(Cell_type, id.outcome) %>%
    summarise(test_count = n(), .groups = "drop")
  if (i %% 500 == 0) message("已读取背景文件：", i, "/", length(mr_files))
}

background_final <- bind_rows(background_list) %>%
  inner_join(disease_type, by = "id.outcome")
result_final <- bind_rows(result_list)
# Every retained positive pair must belong to the eligible test background.
stopifnot(nrow(result_final) == nrow(result))


cell_type_all <- background_final %>%
  group_by(Cell_type) %>% summarise(N = sum(test_count), .groups = "drop")

cell_system_all <- background_final %>%
  group_by(Cell_type, category) %>% summarise(K = sum(test_count), .groups = "drop")

cell_type_pos <- result_final %>%
  group_by(Cell_type) %>% summarise(n = n(), .groups = "drop")

cell_system_pos <- result_final %>%
  group_by(Cell_type, category) %>% summarise(k = n(), .groups = "drop")

# Include cell type-disease system combinations with zero positive pairs.
enrich_result <- expand_grid(Cell_type = sort(unique(na.omit(immune_type$Cell_type))),
                             category = sort(unique(disease_type$category))) %>%
  left_join(cell_type_all, by = "Cell_type") %>%
  left_join(cell_system_all, by = c("Cell_type", "category")) %>%
  left_join(cell_type_pos, by = "Cell_type") %>%
  left_join(cell_system_pos, by = c("Cell_type", "category")) %>%
  mutate(across(c(N, K, n, k), ~replace_na(.x, 0)))
stopifnot(with(enrich_result, all(k <= K & k <= n & n <= N & N - K - n + k >= 0)))

# Use the upper tail of the hypergeometric distribution to test enrichment.
enrich_result <- enrich_result %>%
  rowwise() %>%
  mutate(pvalue = if (N > 0 & K > 0 & K < N) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA_real_,
         ER = if (N > 0 & K > 0 & n > 0) (k / n) / (K / N) else NA_real_) %>% 
  ungroup()

output_dir <- "E:/immune_disease"
dir.create(output_dir, showWarnings = FALSE)
write_csv(enrich_result, file = file.path(output_dir, "enrich_result.csv"))


