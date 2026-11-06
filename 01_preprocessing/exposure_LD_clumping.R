library(TwoSampleMR)
library(ieugwasr)

# Configure paths for the exposure dataset and an ancestry-matched LD reference.
input_dir <- "E:/immune_disease/"
output_dir <- "E:/immune_disease/immune_cells_1e5/"

# Input CSVs must already be instrument-filtered and contain SNP (or rsid) and pval.
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file_path in files) {
  data <- read.csv(file_path)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  colnames(data)[colnames(data) == "SNP"] <- "rsid"
  data_clump <- ld_clump_local(
    dat = data,
    clump_r2 = 0.001,
    clump_kb = 10000,
    clump_p = 1,
    plink_bin = "E:/immune_disease/plink/plink",
    bfile = "E:/immune_disease/data_maf0.01_rs_ref/data_maf0.01_rs_ref"
  )
  colnames(data_clump)[colnames(data_clump) == "rsid"] <- "SNP"
  
  write.table(data_clump,file = paste0(output_dir,file_name,".csv"),row.names = F,col.names = T,sep = ",",quote = F)
}



