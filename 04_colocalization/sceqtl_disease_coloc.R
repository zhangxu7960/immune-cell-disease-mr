args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])  


library(readr)
library(data.table)
library(coloc)
library(susieR)
library(dplyr)
library(ieugwasr)
library(parallel)
library(Matrix)


positive_result <- list.files("/public/home/zhangxu02/immune_disease/sceQTL/14immune_analysis/1e5ivw_0.001_10000kb_positve",full.names = T)



disease_rawdata_path <- "/public/home/zhangxu02/immune_disease/disease/GWAS_all"
sceQTL_rawdata_path <- "/public/home/zhangxu02/immune_disease/sceQTL/14immune_analysis/raw_data"
sceQTL_filterdata_path <- "/public/home/zhangxu02/immune_disease/sceQTL/14immune_analysis/data_fliter_1e5_0.001_10000kb"
output_susie_path <- "/public/home/zhangxu02/immune_disease/colocalization/sceQTL_disease_colocalization/result_susie"
output_abf_path <- "/public/home/zhangxu02/immune_disease/colocalization/sceQTL_disease_colocalization/result_abf"
error_path <- "/public/home/zhangxu02/immune_disease/colocalization/sceQTL_disease_colocalization/error"

bfile <- "/public/home/zhangxu02/immune_disease/LD/1000_EUR/EUR"
plink_bin <- "/public/home/zhangxu02/immune_disease/LD/plink/plink"
bim_snp <- fread("/public/home/zhangxu02/immune_disease/LD/1000_EUR/EUR.bim")






harmonize_immune_disease_data <- function(immune_rawdata_clean,disease_rawdata_clean){
  
  
  
  dat_merge <- inner_join(immune_rawdata_clean,disease_rawdata_clean,by = "SNP",suffix = c(".immune",".disease"))
  
  dat_merge$label_a1 <- ifelse(dat_merge$effect_allele.immune == dat_merge$effect_allele.disease,"yes","no")
  dat_merge$label_a2 <- ifelse(dat_merge$other_allele.immune == dat_merge$other_allele.disease,"yes","no")
  
  dat_merge$label_a1a2 <- ifelse(dat_merge$effect_allele.immune == dat_merge$other_allele.disease,"yes","no")
  dat_merge$label_a2a1 <- ifelse(dat_merge$other_allele.immune == dat_merge$effect_allele.disease,"yes","no")
  
  dat_merge$effect_allele.disease <- ifelse(dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                            dat_merge$effect_allele.immune,
                                            dat_merge$effect_allele.disease)
  dat_merge$other_allele.disease <- ifelse(dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                           dat_merge$other_allele.immune,
                                           dat_merge$other_allele.disease)
  # Flip the disease effect when the alleles are swapped.
  dat_merge$beta.disease <- ifelse(dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                   -dat_merge$beta.disease,
                                   dat_merge$beta.disease)
  dat_merge[,c("label_a1","label_a2","label_a1a2","label_a2a1")] <- NULL
  
  return(dat_merge)
}


prepare_coloc_data <- function(egene_info,disease_rawdata_gwas,disease_id){
  egene_name <- egene_info$exposure
  egene_cell_type <- egene_info$cell_type
  
  egene_data <- read.table(paste0(sceQTL_filterdata_path,"/",egene_cell_type,"_eqtl_table","/",egene_name,".txt"),
                           colClasses = c("A1" = "character","A2" = "character"),# Keep allele columns as character to preserve literal T alleles.
                           header = T,sep = "\t")
  
  
  # Select the SNP with the smallest P value as the lead SNP.
  lead_snp_info <- egene_data %>% dplyr::slice_min(P_VALUE,n = 1,with_ties = FALSE)
  lead_snp_rsid <- lead_snp_info$RSID
  lead_snp_chr <- lead_snp_info$CHR
  lead_snp_pos <- lead_snp_info$POS
  start_pos <- lead_snp_pos - 500000          
  end_pos <- lead_snp_pos + 500000
  
  
  print(paste0("准备共定位数据： ",disease_id,"_",egene_cell_type,"_",egene_name))
  
  
  immune_rawdata_gwas <- fread(paste0(sceQTL_rawdata_path,"/",egene_cell_type,"_eqtl_table.txt.gz")) %>% 
    dplyr::rename(SNP = RSID,chr = CHR,position = POS,effect_allele = A2,other_allele = A1,eaf = A2_FREQ_ONEK1K,pval = P_VALUE) %>% 
    dplyr::select("CELL_ID","CELL_TYPE","SNP","GENE","GENE_ID","chr","position","effect_allele","other_allele","eaf","beta","se","pval") %>% 
    dplyr::mutate(samplesize = 982)        
  
  
  
  # Extract variants within +/-500 kb of the lead SNP.
  immune_rawdata_clean <- immune_rawdata_gwas %>% dplyr::filter(GENE == egene_name, chr == lead_snp_chr, position >= start_pos, position <= end_pos) %>% 
    dplyr::filter(!is.na(beta),!is.na(se),!is.na(eaf),!is.na(SNP)) %>% 
    dplyr::mutate(beta = as.numeric(beta),se = as.numeric(se),eaf = as.numeric(eaf)) %>% 
    dplyr::mutate(z = beta / se) %>% 
    dplyr::group_by(SNP) %>% 
    dplyr::slice_max(order_by = abs(z), n = 1, with_ties = FALSE) %>% 
    dplyr::ungroup()
  disease_rawdata_clean <- disease_rawdata_gwas %>% dplyr::filter(chr == lead_snp_chr, position >= start_pos, position <= end_pos) %>% 
    dplyr::select("SNP","chr","position","effect_allele","other_allele","eaf","beta","se","pval","samplesize") %>%
    dplyr::filter(!is.na(beta),!is.na(se),!is.na(eaf),!is.na(SNP)) %>%                      
    dplyr::mutate(beta = as.numeric(beta),se = as.numeric(se),eaf = as.numeric(eaf)) %>% 
    dplyr::mutate(z = beta / se) %>% 
    dplyr::group_by(SNP) %>%           # Keep the largest absolute Z score for duplicate SNPs.
    dplyr::slice_max(order_by = abs(z), n = 1, with_ties = FALSE) %>% 
    dplyr::ungroup()
  
  
  merged_data <- harmonize_immune_disease_data(immune_rawdata_clean,disease_rawdata_clean)
  
  
  merged_data_2 <- merged_data %>% dplyr::filter(eaf.immune > 0, eaf.immune < 1, eaf.disease > 0, eaf.disease < 1) %>% 
    dplyr::mutate(maf.immune = pmin(eaf.immune, 1 - eaf.immune),
                  maf.disease = pmin(eaf.disease, 1 - eaf.disease))
  
  
  rm(immune_rawdata_gwas)
  gc()
  return(list(data = merged_data_2,lead = lead_snp_info))
}



run_coloc_abf <- function(prepare_data,disease_id,disease_name,output_abf_path){
  if (is.null(prepare_data) || is.null(prepare_data$data) || nrow(prepare_data$data) == 0) return(NULL)
  
  data <- prepare_data$data
  lead <- prepare_data$lead
  egene_name <- lead$GENE
  egene_cell_type <- lead$CELL_ID
  
  
  print(paste0("coloc.abf analysis： ",disease_id,"_",egene_cell_type,"_",egene_name))
  
  
  abf_immune <- list(snp = data$SNP,
                   beta = data$beta.immune,
                   varbeta = data$se.immune^2,
                   N = data$samplesize.immune[1],
                   MAF = data$maf.immune,
                   type = "quant")
  abf_disease <- list(snp = data$SNP,
                     beta = data$beta.disease,
                     varbeta = data$se.disease^2,
                     type = "cc")
  abf_res <- coloc.abf(dataset1 = abf_immune,dataset2 = abf_disease)
  
  
  abf_res_summary <- as.data.frame(t(abf_res[["summary"]])) %>% 
    dplyr::mutate(lead_snp = lead$RSID,
                  chromosome = lead$CHR,
                  pos = lead$POS,
                  Gene = egene_name,
                  cell_type <- egene_cell_type,
                  Disease = disease_name,
                  Disease_id = disease_id)
  
  
  abf_res_results <- abf_res$results %>% dplyr::arrange(desc(SNP.PP.H4)) %>% 
    dplyr::mutate(lead_snp = lead$RSID,
                  chromosome = lead$CHR,
                  pos = lead$POS,
                  Gene = egene_name,
                  cell_type <- egene_cell_type,
                  Disease = disease_name,
                  Disease_id = disease_id)
  if(!dir.exists(paste0(output_abf_path,"/",disease_id)))dir.create(paste0(output_abf_path,"/",disease_id))
  if(!dir.exists(paste0(output_abf_path,"/",disease_id,"/",egene_cell_type)))dir.create(paste0(output_abf_path,"/",disease_id,"/",egene_cell_type))
  write_csv(abf_res_results,file = paste0(output_abf_path,"/",disease_id,"/",egene_cell_type,"/",egene_name,".csv"))
  
  
  return(abf_res_summary)
}


run_coloc_susie <- function(prepare_data,disease_id,disease_name){
  if (is.null(prepare_data) || is.null(prepare_data$data) || nrow(prepare_data$data) == 0) return(NULL)
  
  data <- prepare_data$data
  lead <- prepare_data$lead
  egene_name <- lead$GENE
  egene_cell_type <- lead$CELL_ID
  
  
  print(paste0("coloc.abf analysis： ",disease_id,"_",egene_cell_type,"_",egene_name))
  
  
  if(any(data$SNP %in% bim_snp$V2)){
    LD <- ieugwasr::ld_matrix(
      data$SNP,
      with_alleles = F,
      pop = "EUR",
      bfile = bfile,
      plink_bin = plink_bin
    )
    
    # Match the summary statistics and LD matrix to the same SNP order.
    common <- intersect(data$SNP, colnames(LD))
    data_2 <- data %>%
      filter(SNP %in% common) %>%
      arrange(match(SNP, common))
    LD <- LD[common, common]
    
    
    
    
    
    susie_immune <- susie_rss(z = data_2$z.immune,
                              R = LD,
                              n = data_2$samplesize.immune[1],
                              L = 10)
    susie_disease <- susie_rss(z = data_2$z.disease,
                               R = LD,
                               n = data_2$samplesize.disease[1],
                               L = 10)
    susie_res <- coloc.susie(susie_immune,susie_disease)
    
    
    
    susie_res_summary <- as.data.frame(susie_res$summary)
    if(dim(susie_res_summary)[1] > 0){
      susie_result_summary <- susie_res_summary %>% 
        dplyr::mutate(lead_snp = lead$RSID,
                      chromosome = lead$CHR,
                      pos = lead$POS,
                      Gene = egene_name,
                      cell_type <- egene_cell_type,
                      Disease = disease_name,
                      Disease_id = disease_id)
      
      return(susie_result_summary)
    }else{
      return()
    }
  }else{
    sink(paste0(error_path,"/susie_snp_not_in_bfile.txt"),append = T)
    cat(paste0(disease_id,"_",egene_cell_type,"_",egene_name,"\n"))
    sink()
    return()
  }
}





for (i in task_id) {
  print(paste0("============================",i,"============================="))
  mr_result <- read_csv(positive_result[i])
  disease_id <- mr_result$id.outcome[1]
  disease_name <- mr_result$outcome[1]
  
  
  disease_rawdata_gwas <- fread(paste0(disease_rawdata_path,"/",disease_id,".gz"))
  if(all(c('SNP','CHR','POS','BETA','SE','EA','OA','P','N','EAF') %in% colnames(disease_rawdata_gwas))){
    disease_rawdata_gwas <- disease_rawdata_gwas %>% 
      dplyr::rename(SNP = SNP,
                    chr = CHR,
                    position = POS,
                    beta = BETA,
                    se = SE,
                    pval = P,
                    effect_allele = EA,
                    other_allele = OA,
                    samplesize = N,
                    eaf = EAF)
  }
  required_col <- c("SNP","chr","position","effect_allele","other_allele","eaf","beta","se","pval","samplesize")
  if(!all(required_col %in% colnames(disease_rawdata_gwas))){
    sink(paste0(error_path,"/","disease_col_is_shortage.txt"),append = T)
    cat(disease_id)
    sink()
    next
  }
  
  
  mr_result_list <- split(mr_result,seq_len(nrow(mr_result)))
  immune_disease_data_list <- lapply(mr_result_list,function(x) prepare_coloc_data(x,disease_rawdata_gwas,disease_id))
  
  
  
  

  coloc_susie_result_list <- list()
  coloc_susie_result_list <- lapply(immune_disease_data_list,function(x) run_coloc_susie(x,disease_id,disease_name))

  coloc_susie_result_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, coloc_susie_result_list)

  coloc_susie_result <- do.call(rbind, coloc_susie_result_list)
  if(!is.null(coloc_susie_result) && nrow(coloc_susie_result) > 0){
    write_csv(coloc_susie_result,file = paste0(output_susie_path,"/",disease_id,".csv"))
  }else{
    sink(file = paste0(error_path,"/no_susie_result.txt"),append = T)
    cat(paste0(disease_id," \n"))
    sink()
  }
  
  
  print(paste0("Completed the co-localization analysis： ",disease_id))
  print("==============================================================")
}



