library(readr)
library(grid)
library(ragg)

result_path <- "E:/immune_disease/dementia_monocytes/association_network"
edges <- read_csv(file.path(result_path, "network_associations.csv"), show_col_types = FALSE)
nodes <- read_csv(file.path(result_path, "network_nodes.csv"), show_col_types = FALSE)
stopifnot(nrow(edges) == 19L, nrow(nodes) == 14L)

figure_width_mm <- 95
figure_height_mm <- 110
figure_dpi <- 600
font_family <- "Arial"
negative_color <- "#4575B4"
positive_color <- "#CB8068"
text_color <- "#25313B"

focal_ids <- c("GCST90001991", "GCST90001988", "GCST90002010")
dementia_id <- "DGWAS-105"

immune_order <- c("GCST90001580", "GCST90001988", "GCST90001476",
                  "GCST90001991", "GCST90002010", "GCST90002008", "GCST90001550")
outcome_order <- c("DGWAS-1041", "DGWAS-1565", "DGWAS-2884", "DGWAS-105",
                   "DGWAS-205", "DGWAS-2458", "DGWAS-2727")
stopifnot(setequal(immune_order, nodes$id[nodes$node_type == "Immune phenotype"]),
          setequal(outcome_order, nodes$id[nodes$node_type == "Outcome GWAS"]))

immune_labels <- c(
  GCST90001991 = "HLA-DR on CD14+ monocytes",
  GCST90001988 = "HLA-DR on CD14+ CD16\u2212 monocytes",
  GCST90002010 = "HLA-DR on monocytes",
  GCST90001476 = "HLA-DR++ monocytes % leukocytes",
  GCST90002008 = "CCR2 on monocytes",
  GCST90001580 = "CD14+ CD16+ monocytes AC",
  GCST90001550 = "CM CD8br % T cells")
outcome_labels <- setNames(nodes$raw_label[match(outcome_order, nodes$id)], outcome_order)

wrap_mm <- function(label, width_mm) {
  current <- ""
  lines <- character()
  for (word in strsplit(label, " +")[[1]]) {
    candidate <- if (current == "") word else paste(current, word)
    if (current != "" && convertWidth(stringWidth(candidate), "mm", TRUE) > width_mm) {
      lines <- c(lines, current)
      current <- word
    } else current <- candidate
  }
  paste(c(lines, current), collapse = "\n")
}
invisible(ragg::agg_capture(width = figure_width_mm, height = figure_height_mm,
                             units = "mm", res = 144))

pushViewport(viewport(gp = gpar(fontfamily = font_family, fontsize = 7.2, fontface = "bold")))
immune_wrapped <- vapply(immune_labels, wrap_mm, character(1), width_mm = 28)
outcome_wrapped <- vapply(outcome_labels, wrap_mm, character(1), width_mm = 25.5)
popViewport()
invisible(dev.off())

immune_y <- setNames(seq(16, 88, by = 12), immune_order)
outcome_y <- setNames(seq(16, 88, by = 12), outcome_order)
immune_x <- 34
outcome_x <- 63.5

draw_bipartite_network <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, figure_width_mm),
                        yscale = c(figure_height_mm, 0)))
  u <- function(x) unit(x, "native")
  txt <- function(label, x, y, size = 7.2, just = "left", face = "plain", color = text_color) {
    grid.text(label, u(x), u(y), just = just,
               gp = gpar(fontfamily = font_family, fontsize = size,
                         fontface = face, col = color, lineheight = 1.10))
  }
  rect <- function(x, y, w, h, fill = NA, border = NA, lwd = 0.5) {
    grid.rect(u(x), u(y), width = u(w), height = u(h),
                gp = gpar(fill = fill, col = border, lwd = lwd))
  }
  seg <- function(x0, y0, x1, y1, color, lty = "solid", lwd = 0.8) {
    grid.segments(u(x0), u(y0), u(x1), u(y1),
                    gp = gpar(col = color, lty = lty, lwd = lwd))
  }
  txt("Immune phenotypes", 30.5, 5, 7.6, just = "right", face = "bold")
  txt("Disease outcomes", 67, 5, 7.6, face = "bold")

  
  draw_order <- order(edges$is_dementia)
  t <- seq(0, 1, length.out = 80)
  curve_x <- (1 - t)^3 * immune_x + 3 * (1 - t)^2 * t * (immune_x + 14) +
    3 * (1 - t) * t^2 * (outcome_x - 14) + t^3 * outcome_x
  for (i in draw_order) {
    y0 <- immune_y[edges$id.exposure[i]]
    y1 <- outcome_y[edges$id.outcome[i]]
    curve_y <- ((1 - t)^3 + 3 * (1 - t)^2 * t) * y0 +
      (3 * (1 - t) * t^2 + t^3) * y1
    color <- if (edges$b[i] > 0) positive_color else negative_color
    grid.lines(u(curve_x), u(curve_y),
                 gp = gpar(col = adjustcolor(color, alpha.f = if (edges$is_dementia[i]) 1 else 0.55),
                           lwd = 0.8, lty = "solid"))
  }
  for (id in immune_order) {
    focal <- id %in% focal_ids
    grid.circle(u(immune_x), u(immune_y[id]), r = unit(1.35, "mm"),
                   gp = gpar(fill = if (focal) negative_color else "#D7E1ED",
                             col = if (focal) negative_color else "#9BAAB9", lwd = 0.6))
    txt(immune_wrapped[id], 30.5, immune_y[id], 7.2, "right",
          if (focal) "bold" else "plain")
  }
  for (id in outcome_order) {
    focal <- id == dementia_id
    rect(outcome_x, outcome_y[id], 2.7, 2.7,
           if (focal) "#F8E6D5" else "#E9EEF4",
           if (focal) positive_color else "#9BAAB9", 0.6)
    txt(outcome_wrapped[id], 67, outcome_y[id], 7.2,
          face = if (focal) "bold" else "plain")
  }
  seg(5, 98, 11, 98, negative_color)
  txt("\u03b2 < 0", 14, 98, 6.7)
  seg(55, 98, 61, 98, positive_color)
  txt("\u03b2 > 0", 64, 98, 6.7)
  grid.circle(u(6.5), u(105), r = unit(1.1, "mm"),
                 gp = gpar(fill = negative_color, col = NA))
  txt("Focal monocyte traits", 10, 105, 6.3)
  rect(58, 105, 2.2, 2.2, "#F8E6D5", positive_color)
  txt("Dementia", 62, 105, 6.3)

  popViewport()
}

stem <- file.path(result_path, "Dementia_bipartite_network")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = figure_width_mm / 25.4,
                     height = figure_height_mm / 25.4, family = font_family)
tryCatch(draw_bipartite_network(), finally = dev.off())
ragg::agg_png(paste0(stem, ".png"), width = figure_width_mm, height = figure_height_mm,
              units = "mm", res = figure_dpi, background = "white")
tryCatch(draw_bipartite_network(), finally = dev.off())
ragg::agg_tiff(paste0(stem, ".tif"), width = figure_width_mm, height = figure_height_mm,
               units = "mm", res = figure_dpi, compression = "lzw", background = "white")
tryCatch(draw_bipartite_network(), finally = dev.off())
