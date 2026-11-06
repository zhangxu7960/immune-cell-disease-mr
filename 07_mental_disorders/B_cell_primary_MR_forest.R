library(readr)
library(dplyr)
library(grid)
library(forestploter)
library(ragg)


result_path <- "E:/immune_disease/mental_disorders_B_cells/primary_MR_forest"
figure_specs <- list(
  list(data_name = "B_cell_primary_MR_main_DGWAS-166_2628_6186_481.csv",
       figure_name = "B_cell_primary_MR_forest",
       expected_ids = c("DGWAS-166", "DGWAS-2628", "DGWAS-6186", "DGWAS-481"),
       expected_role = "Main figure"),
  list(data_name = "B_cell_primary_MR_supplement_DGWAS-4147.csv",
       figure_name = "B_cell_primary_MR_forest_supplement_DGWAS-4147",
       expected_ids = "DGWAS-4147",
       expected_role = "Supplementary figure")
)



figure_width_mm <- 170
layout_height_mm <- 85               
dataset_note_height_mm <- 5
figure_dpi <- 600
base_size <- 8
row_height_mm <- 5.5                   
font_family <- "Arial"
point_color <- "#4575B4"
refline_color <- "#777777"
background_colors <- c("#FFFFFF", "#F3F5F6")


# Plot the original ORs and confidence limits and retain the full reported values.
x_limits <- c(0.50, 1.10)
x_ticks <- c(0.50, 0.75, 1.00)



make_primary_forest <- function(data_file, figure_name, expected_ids, expected_role) {
plot_result <- read_csv(data_file, show_col_types = FALSE) %>%
  arrange(plot_order)
required_columns <- c("exposure", "disease_label", "id.exposure", "id.outcome",
                      "b", "se", "or", "or_lci95", "or_uci95", "FDR",
                      "cell_type", "category", "method", "plot_order",
                      "figure_role", "B_cell_eGene_MR_association_n")

dataset_note <- paste0(if (expected_role == "Main figure") "Datasets: " else
                         "Supplementary dataset: ",
                       paste(expected_ids, collapse = ", "))
ci_left_outside <- which(plot_result$or_uci95 < x_limits[1])
ci_right_outside <- which(plot_result$or_lci95 > x_limits[2])
ci_truncated <- which(plot_result$or_lci95 < x_limits[1] |
                       plot_result$or_uci95 > x_limits[2])



outcome_labels <- c("Anxiety disorders" = "Anxiety disorder",
                    "Mental disorder" = "Mental disorder",
                    "Depression-related phenotype" = "Depression-related",
                    "Vascular dementia" = "Vascular dementia")
outcome_display <- unname(outcome_labels[plot_result$disease_label])

fdr_display <- ifelse(plot_result$FDR < 0.001,
                      "<0.001",
                      sprintf("%.3f", plot_result$FDR))

forest_data <- data.frame(
  "Immune phenotype" = plot_result$exposure,
  "Outcome" = outcome_display,
  " " = rep(paste(rep(" ", 18), collapse = ""), nrow(plot_result)),
  "FDR" = fdr_display,
  "OR (95% CI)" = sprintf("%.3f (%.3f to %.3f)",
                           plot_result$or, plot_result$or_lci95,
                           plot_result$or_uci95),
  check.names = FALSE
)



invisible(ragg::agg_capture(width = figure_width_mm, height = layout_height_mm,
                            units = "mm", res = 144, background = "white"))
layout_device <- dev.cur()

forest_style <- forest_theme(
  base_size = base_size,
  base_family = font_family,
  ci_pch = 15,
  ci_col = point_color,
  ci_fill = point_color,
  ci_lwd = 1,
  ci_Theight = 0.12,
  refline_gp = gpar(col = refline_color, lty = "dashed", lwd = 0.8),
  xaxis_gp = gpar(fontsize = base_size - 0.3, fontfamily = font_family,
                 col = "#333333", lwd = 0.6),
  xlab_gp = gpar(fontsize = base_size, fontfamily = font_family),
  core = list(fg_params = list(hjust = 0, x = 0.03),
              bg_params = list(fill = background_colors, col = NA),
              padding = unit(c(1.0, 1.0), "mm")),
  colhead = list(fg_params = list(hjust = 0, x = 0.03, fontface = "bold"),
                 bg_params = list(fill = "#FFFFFF", col = NA),
                 padding = unit(c(2.0, 1.0), "mm"))
)


p <- forest(
  forest_data,
  est = plot_result$or,
  lower = plot_result$or_lci95,
  upper = plot_result$or_uci95,
  sizes = 0.18,
  ci_column = 3,
  ref_line = 1,
  x_trans = "log",
  xlim = x_limits,
  ticks_at = x_ticks,
  ticks_digits = 3,
  xlab = "Odds ratio (log scale)",
  theme = forest_style
)


add_outside_arrow <- function(plot, rows, side) {
  for (i in rows) {
    arrow_x <- if (side == "left") c(0.13, 0.01) else c(0.87, 0.99)
    plot <- add_grob(
      plot, row = i, col = 3, part = "body", order = "top",
      gb_fn = segmentsGrob,
      x0 = unit(arrow_x[1], "npc"), x1 = unit(arrow_x[2], "npc"),
      y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
      arrow = arrow(length = unit(1.3, "mm"), type = "open", ends = "last"),
      gp = gpar(col = point_color, lwd = 1),
      name = paste0("offscale-", side, "-", i)
    )
  }
  plot
}
p <- add_outside_arrow(p, ci_left_outside, "left")
p <- add_outside_arrow(p, ci_right_outside, "right")


p <- edit_plot(p, col = 4, which = "text", x = unit(1.8, "mm"), hjust = 0)
p <- edit_plot(p, col = 4, part = "header", which = "text",
                x = unit(1.8, "mm"), hjust = 0)


group_start <- which(!duplicated(plot_result$disease_label))
for (i in group_start[group_start > 1L]) {
  p <- add_border(p, row = i, part = "body", where = "top",
                  gp = gpar(col = "#D9DDE0", lwd = 0.7))
}
p <- add_border(p, part = "header", where = "bottom",
                gp = gpar(col = "#777777", lwd = 0.7))



column_width_mm <- c(45, 34, 40, 13, 34)
table_columns <- sort(unique(p$layout$l[p$layout$name == "colhead-fg"]))

p$widths[table_columns] <- unit(column_width_mm / sum(column_width_mm) *
                                (figure_width_mm - 4), "mm")
p$widths[c(1, length(p$widths))] <- unit(2, "mm")
p$heights[c(1, length(p$heights))] <- unit(2, "mm")


body_rows <- sort(unique(p$layout$t[p$layout$name == "core-fg"]))
p$heights[body_rows] <- unit.pmax(p$heights[body_rows], unit(row_height_mm, "mm"))

plot_height_mm <- convertHeight(sum(p$heights), "mm", valueOnly = TRUE)

figure_height_mm <- plot_height_mm + dataset_note_height_mm + 1
note_width_mm <- convertWidth(grobWidth(textGrob(dataset_note,
  gp = gpar(fontfamily = font_family, fontsize = 6.5))), "mm", valueOnly = TRUE)
invisible(dev.off(layout_device))



save_forest <- function(plot, output_path = result_path) {
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  draw_forest <- function() {
    grid.newpage()
    grid.text(dataset_note, x = unit(3.5, "mm"),
              y = unit(1, "npc") - unit(2.2, "mm"),
              just = c("left", "center"),
              gp = gpar(fontfamily = font_family, fontsize = 6.5,
                        col = "#5D6772"))
    pushViewport(viewport(x = 0.5,
                          y = unit(1, "npc") - unit(dataset_note_height_mm, "mm"),
                          width = sum(plot$widths), height = sum(plot$heights),
                          just = c("center", "top")))
    grid.draw(plot)
    popViewport()
  }

  pdf_file <- file.path(output_path, paste0(figure_name, ".pdf"))
  grDevices::cairo_pdf(pdf_file, width = figure_width_mm / 25.4,
                       height = figure_height_mm / 25.4, family = font_family)
  tryCatch(draw_forest(), finally = dev.off())

  png_file <- file.path(output_path, paste0(figure_name, ".png"))
  ragg::agg_png(png_file, width = figure_width_mm, height = figure_height_mm,
                units = "mm", res = figure_dpi, background = "white")
  tryCatch(draw_forest(), finally = dev.off())

  invisible(c(pdf = pdf_file, png = png_file))
}


save_forest(p)
invisible(list(plot = p, data = plot_result,
               width_mm = figure_width_mm, height_mm = figure_height_mm,
               dataset_note = dataset_note))
}


figure_outputs <- lapply(figure_specs, function(spec) {
  make_primary_forest(file.path(result_path, spec$data_name),
                      spec$figure_name, spec$expected_ids, spec$expected_role)
})

