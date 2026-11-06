if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"))
}

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
output_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
local_libraries <- file.path(c(Sys.getenv("LOCALAPPDATA"), "C:/Users/23883/AppData/Local"),
  "R", "win-library", "4.3")
.libPaths(unique(c(local_libraries[dir.exists(local_libraries)], .libPaths())))
required <- c("circlize", "ComplexHeatmap", "RColorBrewer", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

data_file <- file.path(output_dir, "plot_data.csv")
if (!file.exists(data_file)) {
  stop("Run the data-preparation script in this folder first.")
}
dat <- read.csv(data_file, check.names = FALSE, fileEncoding = "UTF-8")
dat <- dat[order(dat$plot_order), ]
categories <- unique(dat[c("category", "sector_order")])
categories <- categories[order(categories$sector_order), ]
stopifnot(nrow(dat) == 165L, nrow(categories) == 19L, !anyDuplicated(dat[c("id.exposure", "id.outcome")]))
method_columns <- c("MR_Egger", "Simple_mode", "Simple_median", "IVW")
method_labels <- c("MR Egger", "Simple mode", "Simple median", "IVW")
stopifnot(all(is.finite(as.matrix(dat[method_columns]))))
split <- factor(dat$category, levels = categories$category)
category_matrix <- matrix(dat$category, ncol = 1)
rownames(category_matrix) <- paste(dat$id.exposure, dat$id.outcome, sep = "::")
category_colors <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel1"))(nrow(categories)),
  categories$category)
col_beta <- circlize::colorRamp2(c(-3, -1, 0, 1, 3), c("#2C3E78", "#5CA8D3", "#F5F5F5", "#FFB3A7",
  "#BC3C28"))

# Asterisks identify |Beta| > 1, not P-value significance.
figure_width_mm <- 170
figure_height_mm <- 135
figure_dpi <- 600
font_family <- "Arial"
main_width_fraction <- 0.64
category_track_height <- 0.065
method_track_height <- 0.075
category_legend_fontsize <- 6.5
category_legend_lineheight <- 1
category_legend_row_gap_mm <- 0.25
category_legend_swatch_mm <- 2.2
category_legend_swatch_gap_mm <- 1.6
category_legend_text_gap_mm <- 0.8
category_legend_heading_gap_mm <- 1.3
# Use a generic font family for off-screen legend metrics.
legend_font_family <- "sans"

add_layer_label <- function(label) {
  circos.text(0, 0.5, label, sector.index = levels(split)[1], track.index = get.current.track.index(),
    facing = "bending.inside", niceFacing = TRUE, adj = c(1.15, 0.5), cex = 0.94, font = 2)
}

# Wrap full category names to the available physical text width.
wrap_label_mm <- function(label, max_width_mm, label_gp) {
  words <- strsplit(label, " +")[[1]]
  result <- character()
  current <- ""
  for (word in words) {
    proposed <- if (nzchar(current))
      paste(current, word) else word
    width_mm <- convertWidth(grobWidth(textGrob(proposed, gp = label_gp)), "mm", valueOnly = TRUE)
    if (nzchar(current) && width_mm > max_width_mm) {
      result <- c(result, current)
      current <- word
    } else current <- proposed
  }
  paste(c(result, current), collapse = "\n")
}

draw_category_legend <- function() {
  left_mm <- main_width_fraction * figure_width_mm + 0.2
  swatch_mm <- category_legend_swatch_mm
  text_gap_mm <- category_legend_text_gap_mm
  text_left_mm <- left_mm + swatch_mm + text_gap_mm
  text_width_mm <- figure_width_mm - text_left_mm - 1
  label_gp <- gpar(fontfamily = font_family, fontsize = category_legend_fontsize, lineheight = category_legend_lineheight,
    col = "#222222")
  labels <- vapply(categories$category, wrap_label_mm, character(1), max_width_mm = text_width_mm,
    label_gp = label_gp)
  stopifnot(all(gsub("[[:space:]]+", " ", labels) == categories$category))
  heights <- vapply(labels, function(label) {
    max(swatch_mm, convertHeight(grobHeight(textGrob(label, gp = label_gp)), "mm", valueOnly = TRUE))
  }, numeric(1))
  row_gap_mm <- category_legend_row_gap_mm
  adjacent_heights_mm <- (head(heights, -1L) + tail(heights, -1L))/2
  row_pitch_mm <- swatch_mm + category_legend_swatch_gap_mm
  if (any(row_pitch_mm - adjacent_heights_mm < row_gap_mm - 1e-08)) {
    stop("Legend text does not fit the compact fixed row pitch. Adjust text wrapping first.")
  }
  body_height_mm <- heights[1]/2 + row_pitch_mm * (length(labels) - 1L) + tail(heights, 1)/2
  heading_height_mm <- 3.4
  heading_gap_mm <- category_legend_heading_gap_mm
  total_height_mm <- heading_height_mm + heading_gap_mm + body_height_mm
  stopifnot(total_height_mm <= figure_height_mm - 8)
  top_mm <- (figure_height_mm + total_height_mm)/2
  grid.text("Disease categories", x = unit(left_mm, "mm"), y = unit(top_mm, "mm"), just = c("left",
    "top"), gp = gpar(fontfamily = font_family, fontsize = 8, fontface = "bold"))
  first_center_mm <- top_mm - heading_height_mm - heading_gap_mm - heights[1]/2
  centers_mm <- first_center_mm - (seq_along(labels) - 1L) * row_pitch_mm
  stopifnot(all(abs(diff(centers_mm) + row_pitch_mm) < 1e-08), all(row_pitch_mm - adjacent_heights_mm >=
    row_gap_mm - 1e-08))
  for (i in seq_along(labels)) {
    center_mm <- centers_mm[i]
    grid.rect(x = unit(left_mm + swatch_mm/2, "mm"), y = unit(center_mm, "mm"), width = unit(swatch_mm,
      "mm"), height = unit(swatch_mm, "mm"), gp = gpar(fill = category_colors[i], col = "#AAAAAA",
      lwd = 0.35))
    grid.text(labels[i], x = unit(text_left_mm, "mm"), y = unit(center_mm, "mm"), just = c("left",
      "center"), gp = label_gp)
  }
  cat(sprintf("Legend: equal swatch center spacing %.1f mm; equal edge gaps %.1f mm.\n",
    row_pitch_mm, row_pitch_mm - swatch_mm))
}

draw_figure <- function() {
  circos.clear()
  on.exit(circos.clear(), add = TRUE)
  par(fig = c(0, main_width_fraction, 0, 1), mar = c(0, 0, 0, 0), family = font_family, ps = 8,
    xaxs = "i", yaxs = "i")
  circos.par(gap.after = c(rep(2, nrow(categories) - 1L), 40), start.degree = 100, canvas.xlim = c(-1.12,
    1.12), canvas.ylim = c(-1.12, 1.12), track.margin = c(0.008, 0.008), cell.padding = c(0,
    0, 0, 0), points.overflow.warning = FALSE)
  circos.heatmap(category_matrix, col = category_colors, split = split, rownames.side = NULL,
    track.height = category_track_height, cluster = FALSE, bg.border = NA)
  add_layer_label("Category")

  for (j in seq_along(method_columns)) {
    mat <- matrix(dat[[method_columns[j]]], ncol = 1, dimnames = list(rownames(category_matrix),
      NULL))
    circos.heatmap(mat, col = col_beta, track.height = method_track_height, cluster = FALSE,
      bg.border = "#222222", bg.lwd = 0.55)
    add_layer_label(method_labels[j])
    if (method_columns[j] == "IVW") {
      circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        vals <- dat$IVW[dat$category == CELL_META$sector.index]
        idx <- which(abs(vals) > 1)
        if (length(idx)) {
          circos.text(idx - 0.5, rep(0.5, length(idx)), "*", cex = 0.68, col = "black",
          font = 2, facing = "inside")
        }
      }, bg.border = NA)
    }
  }

  pushViewport(viewport(x = main_width_fraction/2, y = 0.5, width = main_width_fraction,
    height = 1))
  grid.text("731 immune cell phenotypes\n6,640 disease GWAS datasets\n165 retained MR associations",
    x = unit(0.5, "npc"), y = unit(0.53, "npc"), gp = gpar(fontfamily = font_family, fontsize = 8,
      fontface = "bold", lineheight = 1.3))
  lgd <- Legend(title = "Beta / LogOR", col_fun = col_beta, at = c(-3, -1, 0, 1, 3), direction = "horizontal",
    title_position = "topcenter", legend_width = unit(30, "mm"), grid_height = unit(2.5,
      "mm"), title_gp = gpar(fontfamily = legend_font_family, fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontfamily = legend_font_family, fontsize = 6.5))
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.425, "npc"))
  grid.text("* |Beta| > 1", x = unit(0.5, "npc"), y = unit(0.37, "npc"), gp = gpar(fontfamily = font_family,
    fontsize = 7, fontface = "bold"))
  popViewport()
  draw_category_legend()
}

stem <- file.path(output_dir, "circular_heatmap")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = figure_width_mm/25.4, height = figure_height_mm/25.4,
  family = font_family)
draw_figure()
invisible(dev.off())
ragg::agg_png(paste0(stem, ".png"), width = figure_width_mm, height = figure_height_mm, units = "mm",
  res = figure_dpi, background = "white")
draw_figure()
invisible(dev.off())
cat("Saved vector PDF and 600 dpi PNG with a compact 19-category legend.\n")
