# Run the heatmap statistics script before this plotting script.
# Tile values and colours show the number of significant associations.

library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)

base_dir <- "E:/immune_disease/sceQTL_14_immune_cells"
figure_dir <- file.path(base_dir, "eGene_coverage_and_disease_associations")
data_dir <- file.path(figure_dir, "plot_data")
output_dir <- file.path(figure_dir, "figures", "individual_panels")
if(!dir.exists(output_dir))dir.create(output_dir, recursive = TRUE)
plot_data <- fread(file.path(data_dir, "heatmap_data.csv"), data.table = FALSE)

cell_order <- c("bin", "bmem", "cd4nc", "cd4et", "cd4sox4",
                "cd8nc", "cd8et", "cd8s100b", "nk", "nkr",
                "monoc", "mononc", "dc", "plasma")
cell_labels <- c(
  "bin" = "B IN", "bmem" = "B Mem",
  "cd4nc" = "CD4 NC", "cd4et" = "CD4 ET", "cd4sox4" = "CD4 SOX4",
  "cd8nc" = "CD8 NC", "cd8et" = "CD8 ET", "cd8s100b" = "CD8 S100B",
  "nk" = "NK", "nkr" = "NK R", "monoc" = "Mono C", "mononc" = "Mono NC",
  "dc" = "DC", "plasma" = "Plasma"
)

disease_file <- "E:/immune_disease/6640-disease-information.xlsx"
disease_categories <- read_excel(disease_file, sheet = 1) %>%
  pull(category) %>% as.character() %>% trimws() %>% unique()
# Wrap full disease names without abbreviating them.
label_wrap_width <- 50


disease_order <- plot_data %>%
  group_by(disease_category) %>%
  summarise(total = sum(n_associations, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total), disease_category) %>% pull(disease_category)
max_count <- max(plot_data$n_associations, na.rm = TRUE)
plot_data <- plot_data %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_order),
    disease_category = factor(disease_category, levels = rev(disease_order)),
    count_label = ifelse(is.na(n_associations), "\u2013", as.character(n_associations)),
    text_color = ifelse(is.na(n_associations) | n_associations == 0, "#707C86",
                         ifelse(n_associations >= max_count * 0.70, "white", "#263746"))
  )


# Zero indicates tested pairs without significant associations; grey tiles indicate no completed tests.
total_associations <- sum(plot_data$n_associations, na.rm = TRUE)
total_outcomes <- plot_data %>% dplyr::distinct(disease_category, n_outcome_gwas) %>%
  summarise(n = sum(n_outcome_gwas)) %>% pull(n)

p <- ggplot(plot_data, aes(x = cell_type, y = disease_category, fill = n_associations)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = count_label, color = text_color),
            size = 2.25, family = "Arial") +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#DFECF5", "#A5C8E2", "#5C9AC5", "#1F6093", "#103D65"),
    limits = c(0, max_count),
    breaks = pretty(c(0, max_count), n = 4), na.value = "#E4E6E8",
    name = "Significant associations (n)"
  ) +
  scale_x_discrete(labels = cell_labels, expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = label_wrap_width),
                   expand = c(0, 0)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = grid::unit(65, "mm"),
                               barheight = grid::unit(2.5, "mm"))) +
  labs(
    x = NULL, y = NULL,
    title = "Significant eGene\u2013disease associations",
    subtitle = paste0("FDR < 0.05 | ", format(total_associations, big.mark = ","),
                      " associations | ", total_outcomes, " outcome GWAS datasets")
  ) +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "#263746", size = 6.5),
    axis.text.y = element_text(color = "#263746", size = 7, lineheight = 0.95,
                               margin = margin(r = 6)),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.5, color = "#586773", margin = margin(b = 8)),
    legend.position = "bottom",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.margin = margin(t = 3),
    plot.margin = margin(6, 6, 6, 6),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(plot = p, filename = file.path(output_dir, "heatmap.png"),
       width = 230, height = 185, units = "mm", dpi = 600, bg = "white")
ggsave(plot = p, filename = file.path(output_dir, "heatmap.pdf"),
       width = 230, height = 185, units = "mm",
       device = grDevices::cairo_pdf, bg = "white")
