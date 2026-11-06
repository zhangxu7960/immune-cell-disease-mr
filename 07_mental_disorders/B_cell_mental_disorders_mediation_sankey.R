if (.Platform$OS.type == "windows") {
  invisible(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
  user_library <- file.path(Sys.getenv("USERPROFILE"), "AppData/Local/R/win-library",
    paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = "."))
  if (dir.exists(user_library)) .libPaths(c(user_library, .libPaths()))
}
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ragg))
analysis_path <- "E:/immune_disease"
result_path <- file.path(analysis_path, "mental_disorders_B_cells/mediation_sankey")
metabolite_file <- file.path(analysis_path, "mediation_MR/mediation_1400_metabolites/media_result_final.csv")
inflammatory_file <- file.path(analysis_path, "mediation_MR/mediation_91_inflammatory_proteins/inflammatory_final_result.csv")
primary_mr_file <- file.path(analysis_path, "immune_disease_MR/immune_disease_fdr_0.05_ple_het_reverse_20261012.csv")
figure_name <- "B_cell_psychiatric_mediation_sankey"
font_family <- "Arial"
font_pt <- 7.2
figure_width_mm <- 230
figure_dpi <- 600
flow_height_mm <- 2.2
colors <- c(immune = "#6E91B8", metabolite = "#77A494", inflammatory = "#9577B4",
  outcome = "#CB8068", text = "#25313B", secondary = "#64717E", highlight = "#365F89")

read_paths <- function(file, mediator_type, id_col, name_col) {
  dat <- read_csv(file, show_col_types = FALSE)
  needed <- c("id.immune_cell", "immune_cell", "cell_type", id_col, name_col,
    "id.disease", "disease", "category", "b1", "b2", "b", "pval1", "pval2", "pval",
    "mediated_effect", "direct_effect", "mediated_proportion")
  dat %>% mutate(source_row = row_number() + 1L) %>%
    filter(cell_type == "B cell", category == "Mental and behavioural disorders") %>%
    transmute(immune_id = id.immune_cell, immune_phenotype = immune_cell,
      mediator_type = mediator_type, mediator_id = .data[[id_col]],
      mediator = .data[[name_col]], outcome_id = id.disease, outcome = disease,
      beta1 = b1, beta2 = b2, beta = b, mediation_effect = mediated_effect,
      direct_effect = direct_effect, mediated_proportion = mediated_proportion,
      pval1, pval2, pval, source_row)
}
paths <- bind_rows(
  read_paths(metabolite_file, "Metabolite", "id.metabolite", "metabolite"),
  read_paths(inflammatory_file, "Inflammatory protein", "id.inflammatory", "inflammatory"))
numeric_cols <- c("beta1", "beta2", "beta", "mediation_effect", "direct_effect",
                  "mediated_proportion", "pval1", "pval2", "pval")
stopifnot(nrow(paths) == 12L,
  !anyDuplicated(paths[c("immune_id", "mediator_id", "outcome_id")]),
  all(vapply(paths[numeric_cols], function(x) all(is.finite(x)), logical(1))),
  all(paths$pval1 < .01 & paths$pval2 < .01 & paths$pval < .01),
  all(paths$mediation_effect * paths$beta > 0),
  all(paths$mediated_proportion > 0 & paths$mediated_proportion < 1),
  all(abs(paths$beta1 * paths$beta2 - paths$mediation_effect) < 1e-10),
  all(abs(paths$beta - paths$mediation_effect - paths$direct_effect) < 1e-10),
  all(abs(paths$mediation_effect / paths$beta - paths$mediated_proportion) < 1e-10))


outcome_order <- c("DGWAS-6186", "DGWAS-6282", "DGWAS-104", "DGWAS-1605",
    "DGWAS-4823", "DGWAS-7096", "DGWAS-2628", "DGWAS-698", "DGWAS-857",
    "DGWAS-1586", "DGWAS-6128", "DGWAS-4381")
paths <- paths %>% arrange(match(outcome_id, outcome_order), immune_id, mediator_id)


paths$immune_label <- paths$immune_phenotype
paths$mediator_label <- paths$mediator
paths$outcome_label <- paths$outcome

paths$outcome_label[paths$outcome_id == "DGWAS-6186"] <-
  "Depression-related phenotype"
primary <- read_csv(primary_mr_file, show_col_types = FALSE) %>%
  filter(cell_type == "B cell", category == "Mental and behavioural disorders", is.finite(FDR), FDR < .05)
primary_keys <- paste(primary$id.exposure, primary$id.outcome, sep = "|")
paths <- paths %>%
  mutate(retained_primary_MR = paste(immune_id, outcome_id, sep = "|") %in% primary_keys,
    path_id = sprintf("P%02d", row_number()), plot_order = row_number()) %>%
  select(path_id, immune_id, immune_phenotype, immune_label, mediator_type, mediator_id, mediator,
    mediator_label, outcome_id, outcome, outcome_label, beta1, beta2, beta,
    mediation_effect, direct_effect, mediated_proportion, pval1, pval2, pval,
    retained_primary_MR, source_row, plot_order)
stopifnot(sum(paths$mediator_type == "Metabolite") == 11L,
  sum(paths$mediator_type == "Inflammatory protein") == 1L,
  sum(paths$retained_primary_MR) == 1L,
  paths$mediator_id[paths$retained_primary_MR] == "GCST90200673")
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
data_file <- file.path(result_path, "B_cell_psychiatric_mediation_plot_data.csv")
write_excel_csv(paths, data_file, na = "")


immune_left <- 55
immune_right <- 57.8
mediator_left <- 122
mediator_right <- 124.8
outcome_left <- 150
outcome_right <- 152.8
body_top <- 8
row_min_mm <- 7.8
text_gap_mm <- 3.4
wrap_mm <- function(text, width_mm, size_pt = font_pt, bold = FALSE) {
  
  tokenize <- function(x) regmatches(x,
    gregexpr("(?i:immunoglobulin D(?:-positive|-negative)?)|B cells|\\S+", x, perl = TRUE))[[1]]
  words <- tokenize(text)
  width_of <- function(x) convertWidth(grobWidth(textGrob(x,
    gp = gpar(fontfamily = font_family, fontsize = size_pt,
      fontface = if (bold) "bold" else "plain"))), "mm", valueOnly = TRUE)
  result <- character()
  current <- ""
  for (word in words) {
    candidate <- if (current == "") word else paste(current, word)
    measured <- width_of(candidate)
    if (current != "" && measured > width_mm) {
      result <- c(result, current)
      current <- word
    } else current <- candidate
  }
  lines <- c(result, current)
  if (length(lines) > 1L) {
    j <- length(lines)
    previous <- tokenize(lines[j-1L])
    if (length(tokenize(lines[j])) == 1L && length(previous) > 1L) {
      candidate <- paste(tail(previous, 1L), lines[j])
      if (width_of(candidate) <= width_mm) {
        lines[j-1L] <- paste(head(previous, -1L), collapse = " ")
        lines[j] <- candidate
      }
    }
  }
  paste(lines, collapse = "\n")
}

invisible(ragg::agg_capture(width = figure_width_mm, height = 180,
  units = "mm", res = 150, background = "white"))
measure_device <- dev.cur()
paths$immune_text <- vapply(seq_len(nrow(paths)), function(i)
  wrap_mm(paths$immune_label[i], 50.5, bold = paths$retained_primary_MR[i]), character(1))
paths$mediator_text <- vapply(seq_len(nrow(paths)), function(i)
  wrap_mm(paths$mediator_label[i], 59, bold = paths$retained_primary_MR[i]), character(1))
paths$outcome_text <- vapply(seq_len(nrow(paths)), function(i)
  wrap_mm(paths$outcome_label[i], figure_width_mm-outcome_right-5,
    bold = paths$retained_primary_MR[i]), character(1))
text_height_mm <- function(text, bold) {
  convertHeight(grobHeight(textGrob(text,
    gp = gpar(fontfamily = font_family, fontsize = font_pt,
      fontface = if (bold) "bold" else "plain", lineheight = 1.12))), "mm", valueOnly = TRUE)
}
immune_heights <- mapply(text_height_mm, paths$immune_text, paths$retained_primary_MR)
mediator_heights <- mapply(text_height_mm, paths$mediator_text, paths$retained_primary_MR)
outcome_heights <- mapply(text_height_mm, paths$outcome_text, paths$retained_primary_MR)

shared_row_heights <- pmax(row_min_mm, mediator_heights+text_gap_mm, outcome_heights+text_gap_mm)
immune_row_heights <- pmax(row_min_mm, immune_heights+text_gap_mm)
body_height <- max(sum(shared_row_heights), sum(immune_row_heights))
shared_row_heights <- shared_row_heights + (body_height-sum(shared_row_heights))/nrow(paths)
immune_row_heights <- immune_row_heights + (body_height-sum(immune_row_heights))/nrow(paths)
paths$mediator_y <- body_top + cumsum(shared_row_heights) - shared_row_heights/2
paths$outcome_y <- paths$mediator_y
body_bottom <- body_top + body_height

immune_order <- order(tolower(paths$immune_phenotype), paths$immune_id)
paths$immune_y <- NA_real_
ordered_heights <- immune_row_heights[immune_order]
paths$immune_y[immune_order] <- body_top + cumsum(ordered_heights) - ordered_heights/2
normalize_text <- function(x) gsub("[[:space:]]+", " ", trimws(x))
stopifnot(all(normalize_text(paths$immune_text) == normalize_text(paths$immune_label)),
  all(normalize_text(paths$mediator_text) == normalize_text(paths$mediator)),
  all(normalize_text(paths$outcome_text) == normalize_text(paths$outcome_label)),
  all(paths$outcome_label[paths$outcome_id != "DGWAS-6186"] ==
        paths$outcome[paths$outcome_id != "DGWAS-6186"]),
  all(paths$mediator_y-outcome_heights/2 >= body_top),
  all(paths$mediator_y+outcome_heights/2 <= body_bottom),
  n_distinct(paths$immune_id) == nrow(paths), n_distinct(paths$mediator_id) == nrow(paths),
  n_distinct(paths$outcome_id) == nrow(paths))
figure_height_mm <- body_bottom + 1.5
invisible(dev.off(measure_device))

bezier <- function(p0, p1, p2, p3, t) {
  (1-t)^3*p0 + 3*(1-t)^2*t*p1 + 3*(1-t)*t^2*p2 + t^3*p3
}
draw_sankey <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, figure_width_mm), yscale = c(figure_height_mm, 0)))
  grid.rect(gp = gpar(fill = "white", col = NA))
  label <- function(text, x, y, size = font_pt, col = colors[["text"]],
    bold = FALSE, just = "left", halo = FALSE) {
    gp <- gpar(fontfamily = font_family, fontsize = size,
      fontface = if (bold) "bold" else "plain", col = col, lineheight = 1.12)
    if (halo) {
      offsets <- expand.grid(dx = c(-.14, 0, .14), dy = c(-.14, 0, .14))
      for (k in seq_len(nrow(offsets))) {
        grid.text(text, unit(x+offsets$dx[k], "native"), unit(y+offsets$dy[k], "native"),
          just = c(just, "center"), gp = gpar(fontfamily = font_family, fontsize = size,
            fontface = if (bold) "bold" else "plain", col = "white", lineheight = 1.12))
      }
    }
    grid.text(text, unit(x, "native"), unit(y, "native"), just = c(just, "center"), gp = gp)
  }
  ribbon <- function(x0, x1, y0, y1, col, alpha) {
    t <- seq(0, 1, length.out = 100)
    x <- bezier(x0, x0+(x1-x0)*.44, x1-(x1-x0)*.44, x1, t)
    y <- bezier(y0, y0, y1, y1, t)
    grid.polygon(x = unit(c(x, rev(x)), "native"),
      y = unit(c(y-flow_height_mm/2, rev(y+flow_height_mm/2)), "native"),
      gp = gpar(fill = adjustcolor(col, alpha.f = alpha), col = NA))
  }
  node <- function(x0, x1, y, height, col) {
    grid.rect(x = unit((x0+x1)/2, "native"), y = unit(y, "native"),
      width = unit(x1-x0, "native"), height = unit(height, "native"), gp = gpar(fill = col, col = NA))
  }
  label("B-cell phenotypes", 3, 3.4, size = 7.8, bold = TRUE)
  label("Candidate mediators", 65, 3.4, size = 7.8, bold = TRUE)
  label("Outcomes", outcome_right+1.6, 3.4, size = 7.8, bold = TRUE)
  
  for (i in order(paths$retained_primary_MR)) {
    highlight <- paths$retained_primary_MR[i]
    med_col <- colors[[if (paths$mediator_type[i] == "Metabolite") "metabolite" else "inflammatory"]]
    ribbon(immune_right, mediator_left, paths$immune_y[i], paths$mediator_y[i],
      colors[["immune"]], if (highlight) .55 else .22)
    ribbon(mediator_right, outcome_left, paths$mediator_y[i], paths$outcome_y[i],
      med_col, if (highlight) .72 else .33)
  }
  for (i in seq_len(nrow(paths))) {
    highlight <- paths$retained_primary_MR[i]
    med_col <- colors[[if (paths$mediator_type[i] == "Metabolite") "metabolite" else "inflammatory"]]
    node(immune_left, immune_right, paths$immune_y[i], flow_height_mm, colors[["immune"]])
    node(mediator_left, mediator_right, paths$mediator_y[i], flow_height_mm, med_col)
    label(paths$immune_text[i], immune_left-1.5, paths$immune_y[i],
      bold = highlight, just = "right", halo = TRUE,
      col = if (highlight) colors[["highlight"]] else colors[["text"]])
    label(paths$mediator_text[i], mediator_left-1.5, paths$mediator_y[i],
      bold = highlight, just = "right", halo = TRUE,
      col = if (highlight) colors[["highlight"]] else colors[["text"]])
    node(outcome_left, outcome_right, paths$outcome_y[i], flow_height_mm, colors[["outcome"]])
    label(paths$outcome_text[i], outcome_right+1.6, paths$outcome_y[i],
      bold = highlight, col = if (highlight) colors[["highlight"]] else colors[["text"]])
  }
  popViewport()
}
pdf_file <- file.path(result_path, paste0(figure_name, ".pdf"))
grDevices::cairo_pdf(pdf_file, width = figure_width_mm/25.4,
  height = figure_height_mm/25.4, family = font_family)
tryCatch(draw_sankey(), finally = dev.off())
png_file <- file.path(result_path, paste0(figure_name, ".png"))
ragg::agg_png(png_file, width = figure_width_mm, height = figure_height_mm,
  units = "mm", res = figure_dpi, background = "white")
tryCatch(draw_sankey(), finally = dev.off())

