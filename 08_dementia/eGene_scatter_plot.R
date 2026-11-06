library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ragg)
library(grid)


result_path <- "E:/immune_disease/dementia_monocytes/eGene_associations"
data_file <- file.path(result_path, "eGene_plot_data.csv")

outcome_id <- "DGWAS-105"
outcome_titles <- c("DGWAS-105" = "Dementia")
outcome_file_names <- c("DGWAS-105" = "Dementia")

font_family <- "Arial"
figure_dpi <- 600
single_width_mm <- 170
single_height_mm <- 72
single_base_size <- 8
single_gene_size <- 2.45               
fdr_cutoff <- 0.05
random_seed <- 1234

cell_order <- c("bin", "bmem", "cd4et", "cd4nc", "cd4sox4", "cd8et", "cd8nc",
                "cd8s100b", "dc", "monoc", "mononc", "nk", "nkr", "plasma")
cell_labels <- c("B IN", "B MEM", "CD4\nET", "CD4\nNC", "CD4\nSOX4", "CD8\nET",
                 "CD8\nNC", "CD8\nS100B", "DC", "MONO\nC", "MONO\nNC",
                 "NK", "NK R", "Plasma")
point_colors <- c(Negative = "#4575B4", NotSig = "#B8BCC2", Positive = "#CB745C")


plot_result <- read_csv(data_file, show_col_types = FALSE)
required_columns <- c("id.outcome", "exposure", "cell_type", "b", "FDR")

plot_result <- plot_result %>%
  arrange(match(id.outcome, outcome_id), match(cell_type, cell_order), exposure) %>%
  mutate(cell_x = match(cell_type, cell_order),
         status = case_when(FDR >= fdr_cutoff ~ "NotSig",
                            b < 0 ~ "Negative", TRUE ~ "Positive"),
         status = factor(status, levels = names(point_colors)))

set.seed(random_seed)
plot_result$x_plot <- plot_result$cell_x + runif(nrow(plot_result), -0.32, 0.32)


make_egene_plot <- function(id) {
  dat <- plot_result %>% filter(id.outcome == id)
  sig <- dat %>% filter(FDR < fdr_cutoff)
  dat_background <- dat %>% filter(FDR >= fdr_cutoff)
  base_size <- single_base_size
  gene_size <- single_gene_size

  ggplot() +
    annotate("rect", xmin = 9.5, xmax = 11.5, ymin = -Inf, ymax = Inf,
             fill = "#EDF3F8", alpha = 0.8) +
    geom_hline(yintercept = 0, colour = "#AEB5BC", linetype = "dashed",
                linewidth = 0.25) +
    geom_point(data = dat_background, aes(x = x_plot, y = b, colour = status),
                 alpha = 0.15, size = 0.40) +
    geom_point(data = sig, aes(x = x_plot, y = b, colour = status),
                 alpha = 0.95, size = 1.55) +
    geom_text_repel(data = sig, aes(x = x_plot, y = b, label = exposure),
                      family = font_family, fontface = "italic", size = gene_size,
                      colour = "#272727", seed = random_seed,
                      box.padding = 0.25,
                      point.padding = 0.12, min.segment.length = 0,
                      segment.colour = "#868C92", segment.size = 0.20,
                      max.overlaps = Inf, max.time = Inf, max.iter = 50000,
                      force = 1.5,
                      force_pull = 0.4, xlim = c(0.5, 14.5)) +
    scale_colour_manual(values = point_colors, breaks = names(point_colors),
                         drop = FALSE, name = NULL, guide = "none") +
    scale_x_continuous(breaks = seq_along(cell_order),
                        labels = cell_labels,
                        limits = c(0.5, 14.5), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = function(x) pretty(x, n = 3),
                        labels = function(x) sprintf("%.3f", x),
                        expand = expansion(mult = c(0.035, 0.035))) +
    labs(title = unname(outcome_titles[id]), x = "Immune cell type", y = expression(beta)) +
    theme_classic(base_size = base_size, base_family = font_family) +
    theme(plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0,
                                     margin = margin(0, 0, 2, 0)),
           axis.text.x = element_text(size = 6.8,
                                       angle = 0, hjust = 0.5, vjust = 0.5,
                                       colour = "#333333", lineheight = 0.95),
           axis.text.y = element_text(size = 7,
                                       colour = "#333333"),
           axis.title = element_text(size = base_size),
           axis.title.x = element_text(margin = margin(t = 2)),
           axis.title.y = element_text(margin = margin(r = 2)),
           axis.line = element_line(linewidth = 0.25),
           axis.ticks = element_line(linewidth = 0.25),
           axis.ticks.length = unit(0.8, "mm"),
           legend.position = "none",
           plot.margin = margin(2, 2, 1, 1.5))
}



save_plot <- function(plot, filename, width_mm, height_mm) {
  ggsave(file.path(result_path, paste0(filename, ".pdf")), plot = plot,
           width = width_mm, height = height_mm, units = "mm",
           device = grDevices::cairo_pdf, family = font_family, bg = "white")
  ggsave(file.path(result_path, paste0(filename, ".png")), plot = plot,
           width = width_mm, height = height_mm, units = "mm",
           device = ragg::agg_png, dpi = figure_dpi, bg = "white")
}

single_plots <- setNames(lapply(outcome_id, make_egene_plot), outcome_id)
for (id in outcome_id) {
  filename <- paste(id, outcome_file_names[[id]], "eGene_beta", sep = "_")
  save_plot(single_plots[[id]], filename,
              single_width_mm, single_height_mm)
}


