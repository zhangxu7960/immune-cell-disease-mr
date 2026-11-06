library(dplyr)
library(readr)
library(grid)
library(ragg)

result_path <- "E:/immune_disease/dementia_monocytes/functional_enrichment"
enrich_file <- file.path(result_path, "GO_KEGG_significant_enrichment_padjust0.05.tsv")
top_n_per_type <- c(BP = 4L, CC = 4L, MF = 4L, KEGG = 4L)
enrich_fdr_cutoff <- 0.05
font_family <- "Arial"
font_pt <- 7
figure_width_mm <- 200
figure_dpi <- 600
node_colors <- c(eGene = "#9BA8B6", MF = "#9577B4", CC = "#77A494",
                 BP = "#4575B4", KEGG = "#CB8068")

ratio_palette_colors <- c("#FFFFFF", "#DFECF5", "#A5C8E2",
                          "#5C9AC5", "#1F6093", "#103D65")
ratio_palette <- scales::gradient_n_pal(ratio_palette_colors)(seq(0, 1, length.out = 256))

# Select the top N terms within each category by adjusted P value, raw P value and ID.
data_enrich_plot <- read_tsv(enrich_file, show_col_types = FALSE) %>%
  select(ONTOLOGY, ID, Description, GeneRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
  filter(ONTOLOGY %in% names(top_n_per_type),
         is.finite(p.adjust), p.adjust > 0, p.adjust < enrich_fdr_cutoff,
         is.finite(pvalue), pvalue >= 0, pvalue <= 1,
         Count > 0, !is.na(geneID), geneID != "") %>%
  arrange(p.adjust, pvalue, ID) %>%
  group_by(ONTOLOGY) %>% mutate(term_rank = row_number()) %>%
  filter(term_rank <= top_n_per_type[as.character(ONTOLOGY)]) %>%
  ungroup() %>% arrange(p.adjust, pvalue, ID) %>%
  mutate(log10_padjust = -log10(p.adjust),
         GeneRatio_numeric = vapply(strsplit(GeneRatio, "/", fixed = TRUE),
           function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)))
stopifnot(nrow(data_enrich_plot) == 16L, !anyDuplicated(data_enrich_plot$ID),
          identical(as.integer(table(factor(data_enrich_plot$ONTOLOGY,
                                           levels = names(top_n_per_type)))),
                    unname(top_n_per_type)),
          all(diff(data_enrich_plot$p.adjust) >= 0))

gene_sets <- lapply(strsplit(data_enrich_plot$geneID, "/", fixed = TRUE),
                    function(x) unique(trimws(x)))
stopifnot(all(lengths(gene_sets) == data_enrich_plot$Count))
genes <- unique(unlist(gene_sets))
links <- bind_rows(lapply(seq_len(nrow(data_enrich_plot)), function(i) {
  data.frame(gene = gene_sets[[i]], term_index = i)
}))
stopifnot(length(genes) > 0L, nrow(links) == sum(data_enrich_plot$Count),
          !anyDuplicated(links[c("gene", "term_index")]))


gene_left <- 19
gene_right <- 23
term_left <- 150
term_right <- 153
term_text_right <- 148.5
term_text_width <- 72
bubble_left <- term_right
bubble_right <- 184

significance_range <- range(data_enrich_plot$log10_padjust)
significance_padding <- max(diff(significance_range) * 0.08, 0.15)
x_range <- significance_range + c(-significance_padding, significance_padding)
x_ticks <- pretty(x_range, n = 3)
x_ticks <- x_ticks[x_ticks >= x_range[1] & x_ticks <= x_range[2]]
bubble_x <- function(x) bubble_left + (x - x_range[1]) /
  diff(x_range) * (bubble_right - bubble_left)
flow_height_mm <- 1.10
body_top <- 4
gene_gap_mm <- 1.5

wrap_mm <- function(text, max_width_mm) {
  words <- strsplit(text, " +")[[1]]
  lines <- character()
  current <- ""
  for (word in words) {
    candidate <- if (current == "") word else paste(current, word)
    if (current != "" && convertWidth(stringWidth(candidate), "mm",
                                      valueOnly = TRUE) > max_width_mm) {
      lines <- c(lines, current)
      current <- word
    } else current <- candidate
  }
  paste(c(lines, current), collapse = "\n")
}


capture <- ragg::agg_capture(width = figure_width_mm, height = 130,
                             units = "mm", res = 120)
pushViewport(viewport(gp = gpar(fontfamily = font_family, fontsize = font_pt)))

term_labels <- vapply(data_enrich_plot$Description,
                       wrap_mm, character(1), max_width_mm = term_text_width)
popViewport()
invisible(dev.off())
text_lines <- lengths(strsplit(term_labels, "\n"))
node_height <- data_enrich_plot$Count * flow_height_mm
row_height <- pmax(node_height + 1.5, text_lines * 2.55 + 1.0)
row_top <- body_top + c(0, head(cumsum(row_height), -1))
term_y <- row_top + row_height / 2
term_ymin <- term_y - node_height / 2
term_ymax <- term_y + node_height / 2
body_bottom <- body_top + sum(row_height)
figure_height_mm <- body_bottom + 14


gene_degree <- vapply(genes, function(g) sum(links$gene == g), integer(1))
gene_height <- gene_degree * flow_height_mm
gene_stack_height <- sum(gene_height) + (length(genes) - 1) * gene_gap_mm
gene_stack_top <- body_top + (body_bottom - body_top - gene_stack_height) / 2
gene_ymin <- gene_stack_top + c(0, head(cumsum(gene_height + gene_gap_mm), -1))
gene_ymax <- gene_ymin + gene_height
gene_y <- (gene_ymin + gene_ymax) / 2
links <- links %>% mutate(gene_index = match(gene, genes)) %>%
  arrange(gene_index, term_index) %>%
  group_by(gene_index) %>% mutate(gene_slot = row_number()) %>% ungroup() %>%
  group_by(term_index) %>% arrange(gene_index, .by_group = TRUE) %>%
  mutate(term_slot = row_number()) %>% ungroup()
links$start_top <- gene_ymin[links$gene_index] +
  (links$gene_slot - 1) * flow_height_mm
links$end_top <- term_ymin[links$term_index] +
  (links$term_slot - 1) * flow_height_mm
stopifnot(all(links$start_top + flow_height_mm <= gene_ymax[links$gene_index] + 1e-8),
          all(links$end_top + flow_height_mm <= term_ymax[links$term_index] + 1e-8))


ratio_range <- c(0, max(data_enrich_plot$GeneRatio_numeric))
ratio_color <- function(value) {
  idx <- 1 + round(pmin(1, pmax(0, (value - ratio_range[1]) / diff(ratio_range))) * 255)
  ratio_palette[idx]
}
count_diameter <- function(count) 4.2 * sqrt(count / max(data_enrich_plot$Count))

draw_sankey_bubble <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, figure_width_mm),
                        yscale = c(figure_height_mm, 0),
                        gp = gpar(fontfamily = font_family, fontsize = font_pt,
                                  col = "#2D3742")))
  u <- function(x) unit(x, "native")
  text_at <- function(label, x, y, size = font_pt, just = "centre",
                      face = "plain", color = "#2D3742") {
    grid.text(label, x = u(x), y = u(y), just = just,
               gp = gpar(fontfamily = font_family, fontsize = size,
                         fontface = face, col = color, lineheight = 1))
  }
  segment <- function(x0, x1, y0, y1 = y0, color = "#DFE5EC", lwd = 0.45) {
    grid.segments(u(x0), u(y0), u(x1), u(y1), gp = gpar(col = color, lwd = lwd))
  }
  dot <- function(x, y, diameter, color) {
    grid.circle(u(x), u(y), r = unit(diameter / 2, "mm"),
                  gp = gpar(fill = color, col = "white", lwd = 0.25))
  }

  
  t <- seq(0, 1, length.out = 90)
  path_x <- (1 - t)^3 * gene_right +
    3 * (1 - t)^2 * t * (gene_right + 49) +
    3 * (1 - t) * t^2 * (term_left - 49) + t^3 * term_left
  for (k in seq_len(nrow(links))) {
    y0 <- links$start_top[k]
    y1 <- links$end_top[k]
    path_y <- ((1 - t)^3 + 3 * (1 - t)^2 * t) * y0 +
      (3 * (1 - t) * t^2 + t^3) * y1
    grid.polygon(u(c(path_x, rev(path_x))),
                   u(c(path_y, rev(path_y + flow_height_mm))),
                   gp = gpar(fill = adjustcolor("#A8B7C7", alpha.f = 0.27),
                             col = adjustcolor("#A8B7C7", alpha.f = 0.20), lwd = 0.15))
  }

  
  for (tick in x_ticks) {
    segment(bubble_x(tick), bubble_x(tick), body_top, body_bottom,
              color = "#E9EDF2", lwd = 0.4)
  }
  for (y in term_y) segment(bubble_left, bubble_right, y, color = "#F0F3F6", lwd = 0.35)
  grid.rect(x = u(mean(c(bubble_left, bubble_right))),
               y = u(mean(c(body_top, body_bottom))),
               width = u(bubble_right - bubble_left),
               height = u(body_bottom - body_top),
               gp = gpar(fill = NA, col = "#AEB8C4", lwd = 0.5))

  for (g in seq_along(genes)) {
    grid.rect(x = u(mean(c(gene_left, gene_right))), y = u(gene_y[g]),
                 width = u(gene_right - gene_left), height = u(gene_height[g]),
                 gp = gpar(fill = node_colors["eGene"], col = NA))
    text_at(genes[g], gene_left - 1.5, gene_y[g], size = 7,
              just = "right", face = if (genes[g] %in% c("HLA-DQA2", "HSPA1B"))
                "bold.italic" else "italic")
  }
  for (i in seq_len(nrow(data_enrich_plot))) {
    grid.rect(x = u(mean(c(term_left, term_right))), y = u(term_y[i]),
                 width = u(term_right - term_left), height = u(node_height[i]),
                 gp = gpar(fill = node_colors[data_enrich_plot$ONTOLOGY[i]], col = NA))
    
    for (offset in list(c(-0.13, 0), c(0.13, 0), c(0, -0.13), c(0, 0.13))) {
      text_at(term_labels[i], term_text_right + offset[1], term_y[i] + offset[2],
                just = "right", color = "white")
    }
    text_at(term_labels[i], term_text_right, term_y[i], just = "right")
    dot(bubble_x(data_enrich_plot$log10_padjust[i]), term_y[i],
          count_diameter(data_enrich_plot$Count[i]),
          ratio_color(data_enrich_plot$GeneRatio_numeric[i]))
  }

  for (tick in x_ticks) {
    segment(bubble_x(tick), bubble_x(tick), body_bottom, body_bottom + 0.8,
              color = "#647181", lwd = 0.5)
    text_at(sprintf("%.3f", tick), bubble_x(tick), body_bottom + 2.7, size = 6.5)
  }
  text_at(expression(-log[10](italic("Adjust P"))),
            mean(c(bubble_left, bubble_right)), body_bottom + 7.3, size = 7.5)

  
  legend_x <- c(18, 36, 66, 99, 130)
  legend_labels <- c("eGene", "Molecular function", "Cellular component",
                     "Biological process", "Pathways")
  for (j in seq_along(node_colors)) {
    grid.rect(x = u(legend_x[j]), y = u(body_bottom + 9.5),
                 width = unit(2.4, "mm"), height = unit(2.4, "mm"),
                 gp = gpar(fill = node_colors[j], col = NA))
    text_at(legend_labels[j], legend_x[j] + 2, body_bottom + 9.5,
              size = 6.8, just = "left")
  }

  
  body_height <- body_bottom - body_top
  count_y <- body_top + body_height * c(0.22, 0.29, 0.36, 0.43)
  text_at("Count", 188.5, count_y[1] - 5, size = 6.8, just = "left")
  for (j in seq_along(2:5)) {
    count <- (2:5)[j]
    dot(191, count_y[j], count_diameter(count), "#8294AA")
    text_at(as.character(count), 194.3, count_y[j], size = 6.5, just = "left")
  }
  bar_top <- body_top + body_height * 0.63
  bar_bottom <- body_top + body_height * 0.91
  text_at("GeneRatio", 187.5, bar_top - 3, size = 6.5, just = "left")
  bar_y <- seq(bar_top, bar_bottom, length.out = 256)
  for (j in seq_along(bar_y)) {
    grid.rect(x = u(189), y = u(bar_y[j]), width = unit(2.6, "mm"),
                 height = u((bar_bottom - bar_top) / 255 + 0.02),
                 gp = gpar(fill = rev(ratio_palette)[j], col = NA))
  }
  ratio_ticks <- seq(0, floor(ratio_range[2] / 0.05) * 0.05, by = 0.05)
  ratio_ticks <- ratio_ticks[ratio_ticks >= ratio_range[1] & ratio_ticks <= ratio_range[2]]
  for (value in ratio_ticks) {
    y <- bar_bottom - (value - ratio_range[1]) / diff(ratio_range) * (bar_bottom - bar_top)
    segment(190.3, 191.2, y, color = "#8390A0", lwd = 0.45)
    text_at(sprintf("%.3f", value), 192.2, y, size = 6, just = "left")
  }
  popViewport()
}


stem <- file.path(result_path, "enrich_sankey_dotplot_dementia")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = figure_width_mm / 25.4,
                     height = figure_height_mm / 25.4, family = font_family, bg = "white")
draw_sankey_bubble()
invisible(dev.off())
ragg::agg_png(paste0(stem, ".png"), width = figure_width_mm, height = figure_height_mm,
               units = "mm", res = figure_dpi, background = "white")
draw_sankey_bubble()
invisible(dev.off())
ragg::agg_tiff(paste0(stem, ".tif"), width = figure_width_mm, height = figure_height_mm,
                units = "mm", res = figure_dpi, compression = "lzw", background = "white")
draw_sankey_bubble()
invisible(dev.off())
