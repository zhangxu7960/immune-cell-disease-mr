rm(list = ls())
library(readr)
library(forestplot)
library(data.table)
library(dplyr)
outputpath <- "E:/immune_disease/mr_result/1e5ivw_mre_adjP/"
ivw_path <- "E:/immune_disease/mr_result/1e5ivw/"
file <- list.files(ivw_path)
# Apply BH within each outcome across all primary IVW exposure tests.
for (i in file) {
  full_path <- paste0(ivw_path,i)
  data <- read_csv(full_path)
  data <- data[data$method == "Inverse variance weighted", ]
  data$pval_adj_fdr <- p.adjust(data$pval,method = "BH")
  outputpath_file <- paste0(outputpath,i)
  fwrite(data,file = outputpath_file,row.names = F)
}

result <- list()

for (i in file) {
  full_path <- paste0(outputpath,i)
  data <- read_csv(full_path)
  data <- subset(data,pval_adj_fdr < 0.05)
  if(nrow(data) > 0){
    result[[length(result) + 1]] <- data
  }
}
final_result <- bind_rows(result)
fwrite(final_result,file = "E:/immune_disease/mr_result/immune_disease_finalresult_mre.csv",row.names = F)
