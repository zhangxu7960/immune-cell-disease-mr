library(ggplot2)
library(dplyr)
library(readr)

output_dir <- "E:/immune_disease"
enrich_result <- read_csv(file.path(output_dir, "enrich_result.csv"),
                          show_col_types = FALSE)
plot_result <- enrich_result %>% filter(k > 0, is.finite(ER))

category_order <- sort(unique(enrich_result$category))
category_labels <- setNames(
  vapply(category_order, function(x) paste(strwrap(x, width = 62), collapse = "\n"),
         character(1)),
  category_order
)

category_labels["Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified"] <-
  "Symptoms, signs and abnormal clinical and laboratory findings,\nnot elsewhere classified"
category_labels["Injury, poisoning and certain other consequences of external causes"] <-
  "Injury, poisoning and certain other consequences\nof external causes"
category_labels["Factors influencing health status and contact with health services"] <-
  "Factors influencing health status and contact\nwith health services"
category_labels["Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism"] <-
  "Diseases of the blood and blood-forming organs\nand certain disorders involving the immune mechanism"
category_labels["Congenital malformations, deformations and chromosomal abnormalities"] <-
  "Congenital malformations, deformations\nand chromosomal abnormalities"

cell_order <- c("B cell", "Dendritic cell", "Granulocyte", "HSC", "lymphocyte",
                "Monocyte", "Myeloid cell", "Natural killer cell", "NKT cell",
                "T cell", "Treg")
stopifnot(setequal(cell_order, unique(enrich_result$Cell_type)),
          !anyDuplicated(category_labels))

colours_grad <- c("#2C3E78", "#3C6BB3", "#5CA8D3",
                  "#FDDBC7", "#EF8A62", "#B2182B")

p <- ggplot(plot_result, aes(x = Cell_type, y = category)) +
  geom_point(data = subset(plot_result, !sig), aes(size = ER),
             color = "grey50", shape = 16, alpha = 1) +
  geom_point(data = subset(plot_result, sig), aes(size = ER, color = log10_p),
             shape = 16, alpha = 1) +
  scale_color_gradientn(colours = colours_grad) +
  scale_size(range = c(2, 9), name = "ER") +
  scale_x_discrete(limits = cell_order, drop = FALSE) +
  scale_y_discrete(limits = category_order, labels = category_labels, drop = FALSE) +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                    color = "grey25"),
        axis.text.y = element_text(size = 12, lineheight = 0.95,
                                    color = "grey25", margin = margin(r = 5)),
        axis.title = element_text(size = 14),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        plot.margin = margin(8, 10, 8, 8)) +
  labs(x = "Immune Cells Type", y = "Disease Category",
       color = "-log10(pval)", size = "ER")

figure_name <- file.path(output_dir, "enrich_dotplot_mrfdr_0.05_pair")
ggsave(plot = p, filename = paste0(figure_name, ".png"),
       device = ragg::agg_png, width = 13, height = 9.8, units = "in",
       dpi = 300, bg = "white")
ggsave(plot = p, filename = paste0(figure_name, ".pdf"),
       device = cairo_pdf, width = 13, height = 9.8, units = "in",
       bg = "white")
ggsave(plot = p, filename = paste0(figure_name, ".tiff"),
       device = ragg::agg_tiff, width = 13, height = 9.8, units = "in",
       dpi = 600, compression = "lzw", bg = "white")

