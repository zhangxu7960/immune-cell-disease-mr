library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)


analysis_path <- "E:/immune_disease/mediation_MR"
output_path <- file.path(analysis_path, "mediation_distribution")
# Locate the companion statistics script in the same directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else {
  getwd()
}
source(file.path(script_dir, "mediation_distribution_statistics.R"), encoding = "UTF-8")
plot_data <- read.csv(file.path(output_path, "mediation_distribution_plot_data.csv"),
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8")

colour_values <- c("Metabolites" = "#375B7A",
                   "Inflammatory factors" = "#375B7A")
base_family <- "Arial"
base_size <- 7.5


# Allocate row heights according to wrapped label lengths.
make_row_positions <- function(data, wrap_width = 57) {
  rows <- data %>%
    distinct(Item, Order) %>%
    arrange(Order) %>%
    mutate(Label = str_wrap(Item, width = wrap_width),
           Lines = str_count(Label, "\n") + 1,
           Row_height = pmax(1, 0.62 * Lines + 0.25),
           y = sum(Row_height) - cumsum(Row_height) + Row_height / 2)
  list(rows = rows, y_limits = c(-0.05, sum(rows$Row_height) + 0.05))
}

theme_dotplot <- function() {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = "#252B31"),
      axis.text.x = element_text(size = 7, colour = "#555D65"),
      axis.title.x = element_text(size = 7.3, margin = margin(t = 5)),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.line.x = element_line(colour = "#6E7780", linewidth = 0.25),
      axis.ticks.x = element_line(colour = "#6E7780", linewidth = 0.25),
      axis.ticks.length = grid::unit(0.8, "mm"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 8.3, face = "bold", colour = "#252B31",
                                 margin = margin(b = 5)),
      plot.title.position = "plot",
      plot.margin = margin(1.8, 3, 1.8, 3, unit = "mm"),
      legend.position = "none"
    )
}


make_panel <- function(data) {
  row_info <- make_row_positions(data)
  rows <- row_info$rows
  y_limits <- row_info$y_limits
  panel_data <- left_join(data, rows %>% select(Item, y), by = "Item")

  labels_plot <- ggplot(rows, aes(x = 1, y = y, label = Label)) +
    geom_text(hjust = 1, size = base_size / .pt, family = base_family,
              colour = "#252B31", lineheight = 1.03) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_dotplot() +
    theme(axis.text = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.line = element_blank(), axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(1.8, 1, 1.8, 1, unit = "mm"))

  make_count_plot <- function(mediator_type) {
    this_data <- filter(panel_data, Mediator_type == mediator_type)
    this_colour <- unname(colour_values[mediator_type])
    this_fill <- this_colour
    max_count <- max(this_data$Count)
    count_breaks <- pretty(c(0, max_count), n = 4)
    count_breaks <- count_breaks[count_breaks >= 0 & count_breaks <= max_count * 1.1]
    if (length(count_breaks) < 2) count_breaks <- 0:max_count
    axis_max <- max(c(count_breaks, max_count))
    label_offset <- max(axis_max * 0.035, 0.25)
    x_limits <- c(-max(axis_max * 0.025, 0.2), axis_max * 1.18)

    ggplot(this_data, aes(x = Count, y = y)) +
      geom_hline(yintercept = rows$y, colour = "#F0F1F3", linewidth = 0.2) +
      geom_vline(xintercept = count_breaks[count_breaks > 0],
                 colour = "#ECEFF1", linewidth = 0.2) +
      geom_vline(xintercept = 0, colour = "#C2C8CD", linewidth = 0.25) +
      geom_point(shape = 21, size = 1.85, stroke = 0.45,
                 colour = this_colour, fill = this_fill) +
      geom_text(aes(x = Count + label_offset, label = Count), hjust = 0,
                size = 7 / .pt, family = base_family, colour = "#343B42") +
      scale_x_continuous(limits = x_limits, breaks = count_breaks, expand = c(0, 0)) +
      scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      labs(title = mediator_type,
           x = "Candidate pathways (n)", y = NULL) +
      theme_dotplot()
  }

  list(plots = list(labels_plot, make_count_plot("Metabolites"),
                    make_count_plot("Inflammatory factors")),
       height = sum(rows$Row_height) + 2.6)
}

disease_panel <- make_panel(filter(plot_data, Panel == "Disease systems"))
cell_panel <- make_panel(filter(plot_data, Panel == "Immune cell types"))


combined_plot <- wrap_plots(
  c(disease_panel$plots, cell_panel$plots), ncol = 3, byrow = TRUE,
  widths = c(2.15, 1.3, 1.3),
  heights = c(disease_panel$height, cell_panel$height)
) +
  plot_annotation(
    theme = theme(
      text = element_text(family = base_family, colour = "#333333"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(2, 4, 2, 4, unit = "mm")
    )
  )


figure_width_mm <- 210
figure_height_mm <- 195
file_stem <- file.path(output_path, "mediation_distribution_dotplot")

ggsave(paste0(file_stem, ".png"), plot = combined_plot,
       device = ragg::agg_png, width = figure_width_mm, height = figure_height_mm,
       units = "mm", dpi = 600, bg = "white")
ggsave(paste0(file_stem, ".pdf"), plot = combined_plot,
       device = grDevices::cairo_pdf, family = base_family,
       width = figure_width_mm, height = figure_height_mm,
       units = "mm", bg = "white")

message("Exported the mediation pathway distribution plot (PNG and PDF) and its plotting data table.")

