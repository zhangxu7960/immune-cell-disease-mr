if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"))
}

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
output_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
analysis_dir <- dirname(dirname(output_dir))
project_dir <- dirname(analysis_dir)
local_libraries <- file.path(c(Sys.getenv("LOCALAPPDATA"), "C:/Users/23883/AppData/Local"),
  "R", "win-library", "4.3")
.libPaths(unique(c(local_libraries[dir.exists(local_libraries)], .libPaths())))
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the R package readxl first.")

primary_file <- file.path(dirname(output_dir), "immune_disease_fdr_0.05_ple_het_reverse_20261012.csv")
methods_file <- file.path(project_dir, "analysis-20260418", "mr_result", "immune_disease_positive_all_mr_result.csv")
disease_file <- file.path(analysis_dir, "6640-disease-information-20261012.xlsx")
stopifnot(file.exists(primary_file), file.exists(methods_file), file.exists(disease_file))
primary <- read.csv(primary_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
methods <- read.csv(methods_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
diseases <- as.data.frame(readxl::read_xlsx(disease_file))
pair_key <- function(x) paste(x$id.exposure, x$id.outcome, sep = "::")
stopifnot(nrow(primary) == 165L, !anyDuplicated(pair_key(primary)), !anyDuplicated(diseases$GWASID))

# Preserve the complete ICD-10 category names.
category_map <- data.frame(category = c("Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
  "Diseases of the circulatory system", "Congenital malformations, deformations and chromosomal abnormalities",
  "Diseases of the digestive system", "Diseases of the ear and mastoid process", "Endocrine, nutritional and metabolic diseases",
  "Diseases of the eye and adnexa", "Diseases of the genitourinary system", "Factors influencing health status and contact with health services",
  "Certain infectious and parasitic diseases", "Injury, poisoning and certain other consequences of external causes",
  "Mental and behavioural disorders", "Diseases of the musculoskeletal system and connective tissue",
  "Neoplasms", "Diseases of the nervous system", "Pregnancy, childbirth and the puerperium",
  "Diseases of the respiratory system", "Diseases of the skin and subcutaneous tissue", "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified"),
  stringsAsFactors = FALSE)
category_map$sector_order <- seq_len(nrow(category_map))

metadata_index <- match(primary$id.outcome, diseases$GWASID)
stopifnot(!anyNA(metadata_index))
current_category <- diseases$category[metadata_index]
stopifnot(!anyNA(current_category), all(nzchar(current_category)), all(primary$category ==
  current_category), all(current_category %in% category_map$category))
primary$category <- current_category
primary$ICD_code <- diseases$ICD_code[metadata_index]
primary$GWAS_source <- diseases$sources[metadata_index]
map_index <- match(primary$category, category_map$category)
primary$sector_order <- category_map$sector_order[map_index]
primary <- primary[order(primary$sector_order, as.integer(sub("DGWAS-", "", primary$id.outcome)),
  primary$id.exposure), ]
primary$plot_order <- seq_len(nrow(primary))
rownames(primary) <- NULL

method_names <- c("MR Egger", "Simple mode", "Simple median", "Inverse variance weighted")
method_columns <- c("MR_Egger", "Simple_mode", "Simple_median", "IVW")
selected <- methods[pair_key(methods) %in% pair_key(primary) & methods$method %in% method_names,
  ]
method_key <- paste(pair_key(selected), selected$method, sep = "::")
stopifnot(nrow(selected) == nrow(primary) * 4L, !anyDuplicated(method_key))
wide <- primary[c("id.exposure", "id.outcome", "exposure", "outcome", "cell_type", "ICD_code",
  "category", "GWAS_source", "FDR", "sector_order", "plot_order")]
names(wide)[names(wide) == "FDR"] <- "IVW_FDR"
for (j in seq_along(method_names)) {
  one_method <- selected[selected$method == method_names[j], ]
  idx <- match(pair_key(primary), pair_key(one_method))
  stopifnot(!anyNA(idx), all(is.finite(one_method$b[idx])))
  wide[[method_columns[j]]] <- one_method$b[idx]
}
stopifnot(max(abs(wide$IVW - primary$b)) < 1e-12)

stopifnot(length(unique(wide$category)) == 19L, length(unique(wide$id.exposure)) == 122L, length(unique(wide$id.outcome)) ==
  127L)
write.csv(wide, file.path(output_dir, "plot_data.csv"), row.names = FALSE, fileEncoding = "UTF-8",
  na = "")
cat("PASS: 165 unique pairs, 122 immune phenotypes, 127 outcome GWAS datasets,\n", "19 current disease categories, and complete four-method effects.\n",
  "Saved the single core source-data table: plot_data.csv\n")
