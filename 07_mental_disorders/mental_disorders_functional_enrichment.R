library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)

data_path <- "E:/immune_disease/mental_disorders_B_cells"
result_path <- file.path(data_path, "functional_enrichment")
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)


mr_file <- file.path(data_path, "all_cells_significant_eGene_associations.csv")
gene_file <- file.path(data_path, "all_cells_significant_eGenes_unique.csv")
disease_id <- c("DGWAS-166", "DGWAS-2628", "DGWAS-6186",
                "DGWAS-481", "DGWAS-4147")
mr_fdr_cutoff <- 0.05


sceqtl_result <- read_csv(mr_file, show_col_types = FALSE)
egene_result <- read_csv(gene_file, show_col_types = FALSE)


# Use eGenes with FDR < 0.05 across all immune cell types in the five selected outcomes.
sceqtl_result <- sceqtl_result %>%
  dplyr::filter(id.outcome %in% disease_id,
                is.finite(FDR), FDR < mr_fdr_cutoff)


# Count each gene once across cell types and outcomes.
egene_mental_disease <- sceqtl_result %>%
  dplyr::filter(!is.na(exposure), exposure != "") %>%
  dplyr::distinct(exposure) %>% dplyr::pull(exposure) %>% sort()
gene_check <- egene_result %>%
  dplyr::filter(!is.na(eGene), eGene != "") %>%
  dplyr::distinct(eGene) %>% dplyr::pull(eGene)



gene_set <- bitr(egene_mental_disease,
                 fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
gene_set <- gene_set %>% dplyr::distinct(SYMBOL, ENTREZID)


gene_unmapped <- tibble(SYMBOL = egene_mental_disease) %>%
  dplyr::anti_join(gene_set, by = "SYMBOL")
gene_entrez <- unique(gene_set$ENTREZID)


go_mental_disease <- enrichGO(gene = gene_entrez,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             minGSSize = 10,
                             maxGSSize = 500,
                             pool = FALSE,
                             readable = TRUE)

go_mental_disease_result <- as_tibble(go_mental_disease@result)
write_tsv(go_mental_disease_result,
          file = file.path(result_path, "go_mental_disease.tsv"))


kegg_mental_disease <- enrichKEGG(gene = gene_entrez,
                                 organism = "hsa",
                                 keyType = "kegg",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 use_internal_data = FALSE)

kegg_mental_disease <- setReadable(kegg_mental_disease,
                                  OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_mental_disease_result <- as_tibble(kegg_mental_disease@result)
write_tsv(kegg_mental_disease_result,
          file = file.path(result_path, "kegg_mental_disease.tsv"))
