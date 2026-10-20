# Shared publication plotting functions.
# Required package: ggplot2.

pub_palette <- c(
  yellow = "#FBE729",
  climate_red = "#B2182B",
  climate_blue = "#2166AC",
  purple = "#552165",
  indigo = "#4E5F9A",
  teal = "#3A809F",
  coral = "#EB605F",
  light_blue = "#C4D7E5",
  lavender = "#BB9FCB",
  pink = "#E5B7D7",
  pale_purple = "#C8B2D5",
  periwinkle = "#ACAFD3",
  black = "#222222",
  grey = "#777777",
  light_grey = "#D9D9D9",
  very_light_grey = "#F4F4F4"
)

pub_scatter <- list(
  color = "#2166AC",
  ribbon = "#D7E3EE",
  point_size = 1.08,
  point_stroke = 0.20,
  line_width = 0.42
)

# Export dimensions in inches; resolution in dots per inch.
panel_size <- list(
  single_width = 2.25,
  single_height = 1.65,
  figure_width = 5.2,
  figure_height = 7.0,
  dpi = 600
)

pub_theme <- function(base_size = 5.8, base_family = "") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = pub_palette[["black"]]),
      axis.title = ggplot2::element_text(size = base_size + 0.4, face = "plain"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
      axis.text = ggplot2::element_text(size = base_size + 1.2, color = pub_palette[["black"]]),
      axis.line = ggplot2::element_line(linewidth = 0.42, color = pub_palette[["black"]]),
      axis.ticks = ggplot2::element_line(linewidth = 0.36, color = pub_palette[["black"]]),
      axis.ticks.length = grid::unit(1.2, "mm"),
      legend.title = ggplot2::element_text(size = base_size, color = pub_palette[["black"]]),
      legend.text = ggplot2::element_text(size = base_size, color = pub_palette[["black"]]),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size + 0.6, face = "bold"),
      plot.tag = ggplot2::element_text(size = base_size + 2.8, face = "bold"),
      plot.margin = ggplot2::margin(2.5, 3.5, 2.5, 2.5),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

pub_alpha <- function(color, alpha = 0.5) {
  grDevices::adjustcolor(color, alpha.f = alpha)
}

# Fix the plotting area at 50 x 35 mm by default.
fix_pub_panel_size <- function(plot,
                               panel_width = grid::unit(50, "mm"),
                               panel_height = grid::unit(35, "mm")) {
  g <- ggplot2::ggplotGrob(plot)

  panel_cols <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  panel_rows <- unique(g$layout$t[grepl("^panel", g$layout$name)])

  g$widths[panel_cols] <- panel_width
  g$heights[panel_rows] <- panel_height

  g
}

save_pub_plot <- function(plot, filename, width, height, dpi = panel_size$dpi) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

# Save matching PNG and PDF files.
save_pub_pair <- function(plot, stem, width, height, dpi = panel_size$dpi) {
  save_pub_plot(plot, paste0(stem, ".png"), width = width, height = height, dpi = dpi)
  save_pub_plot(plot, paste0(stem, ".pdf"), width = width, height = height, dpi = dpi)
}
