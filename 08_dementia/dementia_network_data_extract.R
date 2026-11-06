library(dplyr)
library(readr)

mr_file <- "E:/immune_disease/immune_disease_MR/immune_disease_fdr_0.05_ple_het_reverse_20261012.csv"
result_path <- "E:/immune_disease/dementia_monocytes/association_network"
dementia_id <- "DGWAS-105"
focal_ids <- c("GCST90001991", "GCST90001988", "GCST90002010")

mr <- read_csv(mr_file, show_col_types = FALSE)
required <- c("id.exposure", "id.outcome", "exposure", "outcome", "cell_type",
               "b", "se", "or", "FDR")
stopifnot(all(required %in% names(mr)),
          !anyDuplicated(mr[c("id.exposure", "id.outcome")]),
          all(is.finite(mr$FDR) & mr$FDR < 0.05))
mr <- mr %>% mutate(source_row = row_number() + 1L,
                     from = paste0("I:", id.exposure),
                     to = paste0("D:", id.outcome))

selected_nodes <- paste0("D:", dementia_id)
repeat {
  incident <- mr$from %in% selected_nodes | mr$to %in% selected_nodes
  expanded <- unique(c(selected_nodes, mr$from[incident], mr$to[incident]))
  if (setequal(expanded, selected_nodes)) break
  selected_nodes <- expanded
}
edges <- mr %>% filter(from %in% selected_nodes, to %in% selected_nodes) %>%
  mutate(or_lci95 = exp(b - 1.96 * se),
         or_uci95 = exp(b + 1.96 * se),
         effect_direction = if_else(b > 0, "Positive", "Negative"),
         is_dementia = id.outcome == dementia_id,
         is_focal_phenotype = id.exposure %in% focal_ids)
stopifnot(nrow(edges) == 19L, length(selected_nodes) == 14L,
          all(is.finite(edges$b) & is.finite(edges$se) & edges$se > 0),
          all(abs(log(edges$or) - edges$b) < 1e-6))

degree <- table(c(edges$from, edges$to))
immune_nodes <- edges %>% distinct(id.exposure, exposure, cell_type) %>%
  transmute(node_id = paste0("I:", id.exposure), id = id.exposure,
             node_type = "Immune phenotype", raw_label = exposure,
             cell_type = cell_type, focal = id.exposure %in% focal_ids)
outcome_nodes <- edges %>% distinct(id.outcome, outcome) %>%
  transmute(node_id = paste0("D:", id.outcome), id = id.outcome,
             node_type = "Outcome GWAS", raw_label = outcome,
             cell_type = NA_character_, focal = id.outcome == dementia_id)
nodes <- bind_rows(immune_nodes, outcome_nodes) %>%
  mutate(degree = as.integer(degree[node_id]))
stopifnot(nrow(immune_nodes) == 7L, nrow(outcome_nodes) == 7L,
          nodes$degree[nodes$id == "GCST90001991"] == 7L,
          nodes$degree[nodes$id == dementia_id] == 3L)


dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
write_excel_csv(edges, file.path(result_path, "network_associations.csv"), na = "")
write_excel_csv(nodes, file.path(result_path, "network_nodes.csv"), na = "")

