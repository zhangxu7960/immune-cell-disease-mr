library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

base_dir <- "E:/immune_disease/sceQTL_14_immune_cells"
figure_dir <- file.path(base_dir, "eGene_coverage_and_disease_associations")
data_dir <- file.path(figure_dir, "plot_data")
output_dir <- file.path(figure_dir, "figures", "individual_panels")
if(!dir.exists(output_dir))dir.create(output_dir, recursive = TRUE)
plot_data <- fread(file.path(data_dir, "eGene_coverage_data.csv"), data.table = FALSE)

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

plot_data <- plot_data %>%
  mutate(cell_type = factor(cell_type, levels = rev(cell_order)),
         instrument_group = factor(instrument_group,
                                    levels = c("Single SNP", "Multiple SNPs")))
cell_totals <- plot_data %>% group_by(cell_type) %>%
  summarise(n_egene = sum(n_egene), .groups = "drop")
total_egene <- sum(cell_totals$n_egene)


# Bar lengths and end labels show eGene counts within each cell type.
p_a <- ggplot(plot_data, aes(x = n_egene, y = cell_type, fill = instrument_group)) +
  geom_col(width = 0.66, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = cell_totals,
    aes(x = n_egene, y = cell_type, label = format(n_egene, big.mark = ",", trim = TRUE)),
    inherit.aes = FALSE, hjust = -0.18, size = 2.25,
    family = "Arial", color = "#263746"
  ) +
  scale_fill_manual(
    values = c("Single SNP" = "#A5C8E2", "Multiple SNPs" = "#1F6093"),
    breaks = c("Single SNP", "Multiple SNPs"),
    labels = c("1 SNP", ">1 SNP"), drop = FALSE,
    name = "Instruments per eGene"
  ) +
  scale_y_discrete(labels = cell_labels, expand = expansion(add = 0.6)) +
  scale_x_continuous(
    limits = c(0, max(cell_totals$n_egene) * 1.24),
    breaks = c(0, 1000, 2000),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
    expand = c(0, 0)
  ) +
  labs(
    x = "Number of eGenes", y = NULL,
    title = "eGene coverage",
    subtitle = paste0(format(total_egene, big.mark = ","), " gene\u2013cell pairs")
  ) +
  guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "#8996A0", linewidth = 0.25),
    axis.ticks.x = element_line(color = "#8996A0", linewidth = 0.25),
    axis.text.x = element_text(color = "#263746", size = 6.5),
    axis.text.y = element_text(color = "#263746", size = 7),
    axis.title.x = element_text(color = "#263746", size = 6.5, margin = margin(t = 5)),
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.5, color = "#586773", margin = margin(b = 8)),
    legend.position = "bottom",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.margin = margin(t = 3),
    plot.margin = margin(6, 6, 6, 6),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(plot = p_a, filename = file.path(output_dir, "eGene_coverage.png"),
       width = 85, height = 185, units = "mm", dpi = 600, bg = "white")
ggsave(plot = p_a, filename = file.path(output_dir, "eGene_coverage.pdf"),
       width = 85, height = 185, units = "mm",
       device = grDevices::cairo_pdf, bg = "white")
