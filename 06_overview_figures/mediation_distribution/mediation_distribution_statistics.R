library(dplyr)
library(tidyr)

analysis_path <- "E:/immune_disease/mediation_MR"
output_path <- file.path(analysis_path, "mediation_distribution")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

metabolite_data <- read.csv(
  file.path(analysis_path, "metabolite_mediation_1400", "media_result_final.csv"),
  stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE
)
inflammatory_data <- read.csv(
  file.path(analysis_path, "inflammatory_mediation_91", "inflammatory_final_result.csv"),
  stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE
)


check_pathways <- function(data, mediator_id) {
  id_columns <- c("id.immune_cell", mediator_id, "id.disease")
  stopifnot(all(c(id_columns, "category", "cell_type") %in% names(data)))
  if (anyNA(data[id_columns]) || any(vapply(data[id_columns], function(x) any(x == ""), logical(1)))) {
    stop("Candidate pathways contain missing IDs. Please check the input results.")
  }
  if (anyDuplicated(data[id_columns])) {
    stop("Duplicate immune phenotype-mediator-outcome GWAS combinations were found. Please check the input results.")
  }
  if (anyNA(data$cell_type) || any(data$cell_type == "")) {
    stop("Unclassified immune cells were found. Please check the input results.")
  }
}

check_pathways(metabolite_data, "id.metabolite")
check_pathways(inflammatory_data, "id.inflammatory")

pathway_data <- bind_rows(
  metabolite_data %>%
    transmute(category, cell_type, Mediator_type = "Metabolites"),
  inflammatory_data %>%
    transmute(category, cell_type, Mediator_type = "Inflammatory factors")
)


# Exclude missing disease categories from disease-system counts; use all records for cell-type counts.
make_count_table <- function(data, group_column, panel_name) {
  valid_data <- data %>%
    filter(!is.na(.data[[group_column]]), .data[[group_column]] != "")

  count_data <- valid_data %>%
    transmute(Item = .data[[group_column]], Mediator_type) %>%
    count(Item, Mediator_type, name = "Count") %>%
    complete(Item, Mediator_type = c("Metabolites", "Inflammatory factors"),
             fill = list(Count = 0L))

  group_order <- count_data %>%
    group_by(Item) %>%
    summarise(Total = sum(Count), .groups = "drop") %>%
    arrange(desc(Total), Item) %>%
    mutate(Order = row_number()) %>%
    select(Item, Order)

  denominators <- data %>%
    group_by(Mediator_type) %>%
    summarise(
      Total_in_panel = sum(!is.na(.data[[group_column]]) & .data[[group_column]] != ""),
      Excluded_unclassified = sum(is.na(.data[[group_column]]) | .data[[group_column]] == ""),
      .groups = "drop"
    )

  count_data %>%
    left_join(group_order, by = "Item") %>%
    left_join(denominators, by = "Mediator_type") %>%
    mutate(Panel = panel_name, .before = 1) %>%
    select(Panel, Item, Mediator_type, Count, Order,
           Total_in_panel, Excluded_unclassified) %>%
    arrange(Order, match(Mediator_type, c("Metabolites", "Inflammatory factors")))
}

disease_plot_data <- make_count_table(pathway_data, "category", "Disease systems")
cell_plot_data <- make_count_table(pathway_data, "cell_type", "Immune cell types")
plot_data <- bind_rows(disease_plot_data, cell_plot_data)


count_check <- plot_data %>%
  group_by(Panel, Mediator_type) %>%
  summarise(Count_sum = sum(Count),
            Total_in_panel = first(Total_in_panel), .groups = "drop")
stopifnot(all(count_check$Count_sum == count_check$Total_in_panel))
stopifnot(!anyNA(plot_data$Count), all(plot_data$Count >= 0))

write.csv(plot_data,
          file.path(output_path, "mediation_distribution_plot_data.csv"),
          row.names = FALSE, fileEncoding = "UTF-8", na = "")
print(count_check)

