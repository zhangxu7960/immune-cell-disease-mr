args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])  


library(purrr)
library(readr)
library(data.table)
library(coloc)
library(susieR)
library(dplyr)
library(ieugwasr)
library(parallel)
library(Matrix)

positive_mr_result <- read.csv("/public/home/zhangxu02/immune_disease/colocalization/immune_disease_gwas_colocalization/immune_disease_fdr_0.05_ple_het_reverse.csv")



immune_rawdata_path <- "/public/home/zhangxu02/immune_disease/731_immune_cells/my_immune_cell_rawdata"
disease_rawdata_path <- "/public/home/zhangxu02/immune_disease/disease/GWAS_all"
immune_filterdata_path <- "/public/home/zhangxu02/immune_disease/731_immune_cells/immune_deleteLD_1e-5"
output_susie_path <- "/public/home/zhangxu02/immune_disease/colocalization/immune_disease_gwas_colocalization/result_susie"
output_abf_path <- "/public/home/zhangxu02/immune_disease/colocalization/immune_disease_gwas_colocalization/result_abf"
error_path <- "/public/home/zhangxu02/immune_disease/colocalization/immune_disease_gwas_colocalization/error"

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



prepare_coloc_data <- function(lead_snp_data,immune_rawdata_gwas,disease_rawdata_gwas,immune_id,disease_id){
  lead_snp_info <- as.data.frame(lead_snp_data)
  
  lead_snp_rsid <- lead_snp_info$SNP
  lead_snp_chr <- lead_snp_info$chr
  lead_snp_pos <- lead_snp_info$pos
  start_pos <- lead_snp_pos - 500000          
  end_pos <- lead_snp_pos + 500000
  
  print(paste0("准备共定位数据： ",immune_id,"_",disease_id,"_",lead_snp_rsid))
  
  # Extract variants within +/-500 kb of the lead SNP.
  immune_rawdata_clean <- immune_rawdata_gwas %>% dplyr::filter(chromosome == lead_snp_chr, pos >= start_pos, pos <= end_pos) %>% 
    dplyr::filter(!is.na(beta),!is.na(se),!is.na(eaf),!is.na(SNP)) %>% 
    dplyr::mutate(beta = as.numeric(beta),se = as.numeric(se),eaf = as.numeric(eaf)) %>% 
    dplyr::mutate(z = beta / se) %>% 
    dplyr::group_by(SNP) %>% 
    dplyr::slice_max(order_by = abs(z), n = 1, with_ties = FALSE) %>% 
    dplyr::ungroup()
  disease_rawdata_clean <- disease_rawdata_gwas %>% dplyr::filter(chromosome == lead_snp_chr, pos >= start_pos, pos <= end_pos) %>% 
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
  
  
  return(list(data = merged_data_2,lead = lead_snp_info))
}



run_coloc_abf <- function(prepare_data,immune_name,immune_id,disease_name,disease_id,output_abf_path){
  if (is.null(prepare_data) || is.null(prepare_data$data) || nrow(prepare_data$data) == 0) return(NULL)
  
  data <- prepare_data$data
  lead <- prepare_data$lead
  
  
  print(paste0("coloc.abf analysis： ",immune_id,"_",disease_id,"_",lead$SNP))
  
  
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
    dplyr::mutate(lead_snp = lead$SNP,
                  chromosome = lead$chr,
                  pos = lead$pos,
                  Immune = immune_name,
                  Immune_id = immune_id,
                  Disease = disease_name,
                  Disease_id = disease_id)
  
  
  abf_res_results <- abf_res$results %>% dplyr::arrange(desc(SNP.PP.H4)) %>% 
    dplyr::mutate(lead_snp = lead$SNP,
                  chromosome = lead$chr,
                  pos = lead$pos,
                  Immune = immune_name,
                  Immune_id = immune_id,
                  Disease = disease_name,
                  Disease_id = disease_id)
  if(!dir.exists(paste0(output_abf_path,"/",immune_id,"_",disease_id)))dir.create(paste0(output_abf_path,"/",immune_id,"_",disease_id))
  write_csv(abf_res_results,file = paste0(output_abf_path,"/",immune_id,"_",disease_id,"/",lead$SNP,".csv"))
  
  
  return(abf_res_summary)
}



run_coloc_susie <- function(prepare_data,immune_name,immune_id,disease_name,disease_id){
  if (is.null(prepare_data) || is.null(prepare_data$data) || nrow(prepare_data$data) == 0) return(NULL)
  
  data <- prepare_data$data
  lead <- prepare_data$lead
  
  
  print(paste0("coloc.susie analysis： ",immune_id,"_",disease_id,"_",lead$SNP))
  
  
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
    
    # Adjust the LD matrix if it has a negative eigenvalue.
    if(min(eigen(LD, symmetric=TRUE)$values) < 0){
      LD <- as.matrix(nearPD(LD, corr = TRUE)$mat)
    }
    
    
    
    
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
        dplyr::mutate(lead_snp = lead$SNP,
                      chromosome = lead$chr,
                      pos = lead$pos,
                      Immune = immune_name,
                      Immune_id = immune_id,
                      Disease = disease_name,
                      Disease_id = disease_id)
      
      return(susie_result_summary)
    }else{
      print(paste0("there is no coloc.susie result： ",immune_id,"_",disease_id,"_",lead$SNP))
      return()
    }
  }else{
    sink(paste0(error_path,"/susie_snp_not_in_bfile.txt"),append = T)
    cat(paste0(immune_id,"_",disease_id,"_",lead$SNP,"\n"))
    sink()
    return()
  }
}





for (i in task_id) {
  print(paste0("============================",i,"============================="))
  
  mr_result_info <- positive_mr_result[i, ]
  immune_id <- mr_result_info$id.exposure
  immune_name <- mr_result_info$exposure
  disease_id <- mr_result_info$id.outcome
  disease_name <- mr_result_info$outcome
  
  
  print(paste0("正在运行：",immune_id,"_",disease_id))
  
  
  immune_rawdata_gwas <- fread(paste0(immune_rawdata_path,"/",immune_id,".tsv.gz")) %>% 
    dplyr::mutate(samplesize = as.numeric(samplesize))
  disease_rawdata_gwas <- fread(paste0(disease_rawdata_path,"/",disease_id,".gz"))
  
  
  required_col <- c("SNP","chr","position","effect_allele","other_allele","eaf","beta","se","pval","samplesize")
  if(!all(required_col %in% colnames(disease_rawdata_gwas))){
    sink(paste0(error_path,"/","disease_col_is_shortage.txt"),append = T)
    cat(disease_id)
    sink()
    next
  }
  disease_rawdata_gwas <- disease_rawdata_gwas %>% dplyr::select("SNP","chr","position","effect_allele","other_allele","eaf","beta","se","pval","samplesize") %>%
    dplyr::rename(chromosome = chr, pos = position)
  
  
  
  lead_snp_infos <- read_csv(paste0(immune_filterdata_path,"/",immune_id,".csv"))
  lead_snp_list <- split(lead_snp_infos, seq_len(nrow(lead_snp_infos)))
  
  
  immune_disease_data_list <- lapply(lead_snp_list,function(x) prepare_coloc_data(x,immune_rawdata_gwas,disease_rawdata_gwas,immune_id,disease_id))
  
  
  
  
  
  
  
  
  
  coloc_susie_result_list <- list()
  coloc_susie_result_list <- lapply(immune_disease_data_list,function(x) run_coloc_susie(x,immune_name,immune_id,disease_name,disease_id))
  
  coloc_susie_result_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, coloc_susie_result_list)
  
  coloc_susie_result <- do.call(rbind, coloc_susie_result_list)
  if(!is.null(coloc_susie_result) && nrow(coloc_susie_result) > 0){
    write_csv(coloc_susie_result,file = paste0(output_susie_path,"/",immune_id,"_",disease_id,".csv"))
  }else{
    print(paste0("there is no coloc.susie result： ",immune_id,"_",disease_id))
    sink(file = paste0(error_path,"/no_susie_result.txt"),append = T)
    cat(paste0(immune_id,"_",disease_id," \n"))
    sink()
  }
  
  
  
  
  
  print(paste0("Completed the co-localization analysis： ",immune_id,"_",disease_id))
  print("==============================================================")
  
}



