library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)


data_path <- "E:/immune_disease/dementia_monocytes"
result_path <- file.path(data_path, "functional_enrichment")
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

mr_file <- file.path(data_path, "all_cells_significant_eGene_associations.csv")
gene_file <- file.path(data_path, "all_cells_significant_eGenes_unique.csv")
disease_id <- c("DGWAS-105")
mr_fdr_cutoff <- 0.05


sceqtl_result <- read_csv(mr_file, show_col_types = FALSE)
egene_result <- read_csv(gene_file, show_col_types = FALSE)


# Use dementia-associated eGenes with FDR < 0.05 across all 14 immune cell types.
sceqtl_result <- sceqtl_result %>%
  dplyr::filter(id.outcome %in% disease_id,
                is.finite(FDR), FDR < mr_fdr_cutoff)

# Count each gene once across cell types and outcomes.
egene_dementia <- sceqtl_result %>%
  dplyr::filter(!is.na(exposure), exposure != "") %>%
  dplyr::distinct(exposure) %>% dplyr::pull(exposure) %>% sort()
gene_check <- egene_result %>%
  dplyr::filter(!is.na(eGene), eGene != "") %>%
  dplyr::distinct(eGene) %>% dplyr::pull(eGene)



gene_set <- bitr(egene_dementia,
                 fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
gene_set <- gene_set %>% dplyr::distinct(SYMBOL, ENTREZID)



gene_unmapped <- tibble(SYMBOL = egene_dementia) %>%
  dplyr::anti_join(gene_set, by = "SYMBOL")
gene_entrez <- unique(gene_set$ENTREZID)


go_dementia <- enrichGO(gene = gene_entrez,
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
go_dementia_result <- as_tibble(go_dementia@result)
write_tsv(go_dementia_result,
          file = file.path(result_path, "go_dementia.tsv"))

options(url.method = "libcurl", timeout = max(180, getOption("timeout")))
kegg_dementia <- enrichKEGG(gene = gene_entrez,
                                 organism = "hsa",
                                 keyType = "kegg",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 use_internal_data = FALSE)
kegg_dementia <- setReadable(kegg_dementia,
                                  OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_dementia_result <- as_tibble(kegg_dementia@result)
write_tsv(kegg_dementia_result,
          file = file.path(result_path, "kegg_dementia.tsv"))


enrich_fdr_cutoff <- 0.05
go_significant <- go_dementia_result %>%
  dplyr::filter(is.finite(p.adjust), p.adjust < enrich_fdr_cutoff) %>%
  dplyr::mutate(Database = "GO")
kegg_significant <- kegg_dementia_result %>%
  dplyr::filter(is.finite(p.adjust), p.adjust < enrich_fdr_cutoff) %>%
  dplyr::mutate(Database = "KEGG", ONTOLOGY = "KEGG")

significant_result <- bind_rows(go_significant, kegg_significant) %>%
  dplyr::select(Database, ONTOLOGY, everything()) %>%
  arrange(p.adjust, pvalue, ID)
write_tsv(significant_result,
          file = file.path(result_path, "GO_KEGG_significant_enrichment_padjust0.05.tsv"))

