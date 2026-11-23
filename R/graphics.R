# Shared internal graphics helpers.

# Apply the package theme and set the title, subtitle and caption of a plot.
# Only the supplied elements are passed to ggplot2::labs(), so a title set
# beforehand survives when `title` is NULL. Extra arguments are forwarded to
# ggplot2::labs() as well (axis and legend titles).
.finish_plot <- function(p, title = NULL, subtitle = NULL, caption = NULL, width = 85, ...) {
  args <- list(...)
  if (!is.null(title)) args$title <- title
  if (!is.null(subtitle)) args$subtitle <- subtitle
  if (!is.null(caption)) args$caption <- .wrap_caption(caption, width)
  out <- p + theme_allostatic()
  if (length(args)) out <- out + do.call(ggplot2::labs, args)
  out
}

# Text colour of a heatmap label: white above the threshold, otherwise primary.
.heatmap_text_colour <- function(values, threshold) {
  ifelse(!is.na(values) & values > threshold, "white", .al_palette$primary)
}

# Tile layer and centred label of a heatmap. The fill scale stays with each
# graphic; only the geometry and the label colour are shared. `label` is the
# name of the column holding the label text and the threshold applies to the
# fill column. Returns a list of layers: the tile and, when a label is given,
# the text.
.heatmap_layer <- function(data, x, y, fill, label = NULL, threshold = NULL,
                           width = NULL, height = NULL, label_size = 3, ...) {
  mapping <- ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])
  tile_args <- list(data = data, mapping = mapping)
  if (!is.null(width)) tile_args$width <- width
  if (!is.null(height)) tile_args$height <- height
  layers <- list(do.call(ggplot2::geom_tile, c(tile_args, list(...))))
  if (!is.null(label)) {
    text_args <- list(data = data, mapping = ggplot2::aes(label = .data[[label]]), size = label_size,
      colour = if (is.null(threshold)) .al_palette$primary else .heatmap_text_colour(data[[fill]], threshold))
    layers[[length(layers) + 1L]] <- do.call(ggplot2::geom_text, text_args)
  }
  layers
}
