# Use identical cell coordinates for the bar plot and heatmap.

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(readxl)

base_dir <- "E:/immune_disease/sceQTL_14_immune_cells"
figure_dir <- file.path(base_dir, "eGene_coverage_and_disease_associations")
data_dir <- file.path(figure_dir, "plot_data")
bar_file <- file.path(data_dir, "eGene_coverage_data.csv")
heatmap_file <- file.path(data_dir, "heatmap_data.csv")
output_dir <- file.path(figure_dir, "figures")
if(!dir.exists(output_dir))dir.create(output_dir, recursive = TRUE)

# The default bar height is the eGene count, grouped by the number of SNP instruments.
# Set bar_metric to "n_snp" to display the SNP count instead.
# n_snp counts cell-gene-SNP combinations, not globally unique SNPs.
bar_metric <- "n_egene"

bar_data <- fread(bar_file, data.table = FALSE)
heatmap_data <- fread(heatmap_file, data.table = FALSE)

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
label_wrap_width <- 50


cell_axis <- function(){
  scale_x_continuous(
    breaks = seq_along(cell_order),
    labels = unname(cell_labels[cell_order]),
    limits = c(0.5, length(cell_order) + 0.5),
    expand = c(0, 0)
  )
}

common_theme <- theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "#263746"),
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.5, color = "#586773", margin = margin(b = 7)),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    plot.background = element_rect(fill = "white", color = NA)
  )


bar_data <- bar_data %>%
  mutate(cell_index = match(cell_type, cell_order),
         instrument_group = factor(instrument_group,
                                    levels = c("Single SNP", "Multiple SNPs")),
         bar_value = .data[[bar_metric]])
bar_totals <- bar_data %>% group_by(cell_type, cell_index) %>%
  summarise(total = sum(bar_value), .groups = "drop")

if(bar_metric == "n_egene"){
  bar_title <- "eGene coverage (n)"
  bar_subtitle <- paste0(format(sum(bar_totals$total), big.mark = ","), " gene\u2013cell pairs")
}else{
  bar_title <- "SNP instrument counts (n)"
  bar_subtitle <- paste0(format(sum(bar_totals$total), big.mark = ","),
                         " gene\u2013cell\u2013SNP combinations")
}

p_bar <- ggplot(bar_data, aes(x = cell_index, y = bar_value, fill = instrument_group)) +
  geom_col(width = 0.68, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = bar_totals,
    aes(x = cell_index, y = total, label = format(total, big.mark = ",", trim = TRUE)),
    inherit.aes = FALSE, vjust = -0.45, size = 2.25,
    family = "Arial", color = "#263746"
  ) +
  cell_axis() +
  scale_y_continuous(
    limits = c(0, max(bar_totals$total) * 1.20),
    breaks = pretty(c(0, max(bar_totals$total)), n = 3),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c("Single SNP" = "#A5C8E2", "Multiple SNPs" = "#1F6093"),
    breaks = c("Single SNP", "Multiple SNPs"), labels = c("1 SNP", ">1 SNP"),
    drop = FALSE, name = "Instruments per eGene"
  ) +
  guides(fill = guide_legend(ncol = 1, title.position = "top")) +
  labs(x = NULL, y = NULL) +
  common_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.5),
    axis.ticks.y = element_line(color = "#8996A0", linewidth = 0.25),
    axis.line.y = element_line(color = "#8996A0", linewidth = 0.25),
    legend.position = "right",
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.margin = margin(t = 0, b = 1),
    plot.margin = margin(6, 6, 2, 6)
  )


disease_order <- heatmap_data %>% group_by(disease_category) %>%
  summarise(total = sum(n_associations, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total), disease_category) %>% pull(disease_category)
max_count <- max(heatmap_data$n_associations, na.rm = TRUE)
heatmap_data <- heatmap_data %>%
  mutate(
    cell_index = match(cell_type, cell_order),
    disease_category = factor(disease_category, levels = rev(disease_order)),
    count_label = ifelse(is.na(n_associations), "\u2013", as.character(n_associations)),
    text_color = ifelse(is.na(n_associations) | n_associations == 0, "#707C86",
                         ifelse(n_associations >= max_count * 0.70, "white", "#263746"))
  )
total_associations <- sum(heatmap_data$n_associations, na.rm = TRUE)
total_outcomes <- heatmap_data %>% dplyr::distinct(disease_category, n_outcome_gwas) %>%
  summarise(n = sum(n_outcome_gwas)) %>% pull(n)
heatmap_title <- "Significant eGene\u2013disease associations"
heatmap_subtitle <- paste0("FDR < 0.05\n", format(total_associations, big.mark = ","),
                            " associations\n", total_outcomes, " outcome GWAS datasets")

p_heatmap <- ggplot(heatmap_data, aes(x = cell_index, y = disease_category, fill = n_associations)) +
  geom_tile(width = 1, color = "white", linewidth = 0.35) +
  geom_text(aes(label = count_label, color = text_color),
            size = 2.25, family = "Arial") +
  scale_color_identity() +
  cell_axis() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = label_wrap_width),
                   expand = c(0, 0)) +
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#DFECF5", "#A5C8E2", "#5C9AC5", "#1F6093", "#103D65"),
    limits = c(0, max_count), breaks = pretty(c(0, max_count), n = 4),
    na.value = "#E4E6E8", name = "Significant\nassociations (n)"
  ) +
  guides(fill = guide_colorbar(direction = "vertical", title.position = "top", title.hjust = 0,
                               barwidth = grid::unit(2.7, "mm"),
                               barheight = grid::unit(42, "mm"))) +
  labs(x = "Immune cell type", y = NULL) +
  common_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y = element_text(size = 7, lineheight = 0.95, margin = margin(r = 6)),
    axis.title.x = element_text(size = 7, margin = margin(t = 5)),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.margin = margin(t = 3),
    plot.margin = margin(3, 6, 6, 6)
  )


# Measure legends on a temporary Cairo device to match PDF font metrics.
get_right_legend <- function(plot){
  device_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(device_file, family = "Arial")
  on.exit({grDevices::dev.off(); unlink(device_file)})
  plot_grob <- ggplotGrob(plot + theme(legend.margin = margin(0, 0, 0, 0),
                                      legend.box.margin = margin(0, 0, 0, 0)))
  legend_index <- match("guide-box-right", plot_grob$layout$name)
  if(is.na(legend_index))stop("没有找到右侧图例。")
  plot_grob$grobs[[legend_index]]
}

make_side_notes <- function(title, subtitle, legend){
  # Preserve separate lines for FDR, association totals and GWAS counts.
  subtitle_lines <- strsplit(subtitle, "\n", fixed = TRUE)[[1]]
  subtitle_text <- paste(vapply(subtitle_lines, str_wrap, character(1), width = 25),
                         collapse = "\n")
  title_grob <- grid::textGrob(
    str_wrap(title, width = 22), x = 0, y = 1, just = c("left", "top"),
    gp = grid::gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold",
                    col = "#263746", lineheight = 1.05)
  )
  subtitle_grob <- grid::textGrob(
    subtitle_text, x = 0,
    y = grid::unit(1, "npc") - grid::grobHeight(title_grob) - grid::unit(2, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontfamily = "Arial", fontsize = 6.5,
                    col = "#586773", lineheight = 1.15)
  )
  block_height <- grid::grobHeight(title_grob) + grid::unit(2, "mm") +
    grid::grobHeight(subtitle_grob) + grid::unit(5, "mm") + grid::grobHeight(legend)
  legend_y <- grid::unit(1, "npc") - grid::grobHeight(title_grob) -
    grid::unit(2, "mm") - grid::grobHeight(subtitle_grob) - grid::unit(5, "mm")
  legend_grob <- grid::editGrob(legend, vp = grid::viewport(
    x = 0, y = legend_y, width = grid::grobWidth(legend),
    height = grid::grobHeight(legend), just = c("left", "top")
  ))
  notes_grob <- grid::grobTree(title_grob, subtitle_grob, legend_grob,
    vp = grid::viewport(x = grid::unit(2, "mm"), y = 0.5,
                         width = grid::unit(1, "npc") - grid::unit(4, "mm"),
                         height = block_height, just = c("left", "center")))
  wrap_elements(full = notes_grob, clip = FALSE)
}

bar_legend <- get_right_legend(p_bar)
heatmap_legend <- get_right_legend(p_heatmap)
bar_notes <- make_side_notes(bar_title, bar_subtitle, bar_legend)
heatmap_notes <- make_side_notes(heatmap_title, heatmap_subtitle, heatmap_legend)
p_bar <- p_bar + theme(legend.position = "none")
p_heatmap <- p_heatmap + theme(legend.position = "none")

p_combined <- p_bar + bar_notes + p_heatmap + heatmap_notes +
  plot_layout(design = "AB\nCD", widths = c(1, 0.26),
              heights = c(1, 3.2), guides = "keep")

ggsave(plot = p_combined, filename = file.path(output_dir, "_heatmap.png"),
       width = 250, height = 195, units = "mm", dpi = 600, bg = "white")
ggsave(plot = p_combined, filename = file.path(output_dir, "heatmap.pdf"),
       width = 250, height = 195, units = "mm",
       device = grDevices::cairo_pdf, bg = "white")
