# Count distinct RSIDs within each cell-gene pair to classify single- and multiple-SNP eGenes.

library(data.table)
library(dplyr)

base_dir <- "E:/immune_disease/sceQTL_14_immune_cells"
input_dir <- file.path(base_dir, "14_immune_cells_preprocessed_1e5_0.001_10000kb_F")
output_dir <- file.path(base_dir, "eGene_coverage_and_disease_associations", "plot_data")
if(!dir.exists(output_dir))dir.create(output_dir, recursive = TRUE)

cell_order <- c("bin", "bmem", "cd4nc", "cd4et", "cd4sox4",
                "cd8nc", "cd8et", "cd8s100b", "nk", "nkr",
                "monoc", "mononc", "dc", "plasma")
cell_dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
cell_names <- sub("_eqtl_table$", "", basename(cell_dirs))
if(!setequal(cell_names, cell_order) || anyDuplicated(cell_names)){
  stop("The type of cells in the directory does not match the preset 14 types of cells. Please check the input.")
}


res_list <- lapply(cell_dirs, function(cell_path){
  cell_type <- sub("_eqtl_table$", "", basename(cell_path))
  files <- list.files(cell_path, pattern = "\\.txt$", full.names = TRUE)
  if(length(files) == 0)stop("The cell folder does not contain a gene txt file.", cell_path)

  n_iv <- vapply(files, function(f){
    dat <- fread(f, select = "RSID", colClasses = "character", showProgress = FALSE)
    if(nrow(dat) == 0 || any(is.na(dat$RSID) | trimws(dat$RSID) %in% c("", ".", "NA"))){
      stop("The genetic file is empty or there are missing RSIDs：", f)
    }
    # Count repeated RSIDs across rounds once without selecting or modifying effect estimates.
    as.integer(n_distinct(trimws(dat$RSID)))
  }, integer(1))

  data.frame(
    cell_type = cell_type,
    instrument_group = ifelse(n_iv == 1, "Single SNP", "Multiple SNPs"),
    n_snp = n_iv, stringsAsFactors = FALSE
  )
})

# Count eGenes as cell-gene pairs and sum distinct SNPs within those pairs.
gene_stat <- bind_rows(res_list)
group_stat <- gene_stat %>%
  group_by(cell_type, instrument_group) %>%
  summarise(n_egene = n(), n_snp = sum(n_snp), .groups = "drop")

plot_data <- expand.grid(
  cell_type = cell_order,
  instrument_group = c("Single SNP", "Multiple SNPs"),
  stringsAsFactors = FALSE
) %>%
  left_join(group_stat, by = c("cell_type", "instrument_group")) %>%
  mutate(n_egene = ifelse(is.na(n_egene), 0L, n_egene),
         n_snp = ifelse(is.na(n_snp), 0L, n_snp))

fwrite(plot_data, file.path(output_dir, "eGene_coverage_data.csv"))

