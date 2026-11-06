library(data.table)
library(foreach)
library(TwoSampleMR)
library(ieugwasr)

# Configure data paths, PLINK and an ancestry-matched LD reference panel.
errorpath <- "E:/immune_disease/colocalization/GWAS_sceQTL_colocalization/error"
inputdata <- "E:/immune_disease/sceQTL/pan_baidu/14_immune_cells/sceQTL14_analysis/1_step1/scQTL"
outputdata <- "E:/immune_disease/colocalization/GWAS_sceQTL_colocalization/sceqtl_1e5_r2_0.001_10000kb"
plinkpath <- "E:/immune_disease/immune_cell_rawdata/plink/plink.exe"
bfilepath <- "E:/immune_disease/immune_cell_rawdata/1000_EUR/EUR"

cell_id <- read.table(file = "H:/sceQTL/sceQTL_analysis/1_step1/cell_id.txt",header = T)
cellid <- as.vector(cell_id$id)

bim_snp <- fread("E:/immune_disease/immune_cell_rawdata/1000_EUR/EUR.bim")
for (i in cellid) {
  print(paste0("analysis:",i))
  eqtldata <- fread(paste0(inputdata,"/",i,".txt.gz"))
  if(!dir.exists(paste0(outputdata,"/",i)))dir.create(paste0(outputdata,"/",i))

  eqtl_filter <- eqtldata[eqtldata$P_VALUE < 1e-5, ]
  gene <- unique(eqtl_filter$GENE)
  for (j in gene) {
    if(!file.exists(paste0(outputdata,"/",i,"/",j,".txt"))){
      geneqtl <- subset(eqtl_filter,eqtl_filter$GENE == j)
      if(nrow(geneqtl) > 1){
        if(any(geneqtl$RSID %in% bim_snp$V2)){
          geneqtl2 <- geneqtl[,c("RSID","P_VALUE")]
          colnames(geneqtl2) <- c("rsid", "pval")

          # LD clumping is performed separately for each gene within each cell type.
          clumdf <- ld_clump_local(dat = geneqtl2, clump_kb = 10000, clump_r2 = 0.001,clump_p=1,
                                   bfile = bfilepath,
                                   plink_bin = plinkpath)

          geneqtl3 <- geneqtl[which(geneqtl$RSID%in%clumdf$rsid),]
          write.table(geneqtl3,file = paste0(outputdata,"/",i,"/",j,".txt"),row.names = F,sep = "\t",quote = F)
          print(paste0(i,"_",j," have completed"))
        }else{
          sink(file = paste0(errorpath,"/no_snp_in_bfile.txt"),append = TRUE)
          cat(paste0(i, "_", j, "\n"))
          sink()
        }
      }else{
        write.table(geneqtl,file = paste0(outputdata,"/",i,"/",j,".txt"),row.names = F,sep = "\t",quote = F)
        print(paste0(i,"_",j," have completed"))
      }
    }else{
      print(paste0(i,"_",j," have already existed"))
    }
  }
}
