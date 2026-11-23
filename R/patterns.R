.risk_matrix <- function(result) matrix(result$components$risk, nrow = nrow(result$scores))

# Effective number of independent indicators (Li and Ji 2005).
.meff <- function(eigenvalues) {
  x <- abs(eigenvalues)
  sum((x >= 1) + (x - floor(x)))
}

#' Risk patterns, co-occurrence and panel redundancy
#'
#' A count score hides which biomarkers are at risk: a load of 2 may be two
#' blood-pressure readings (one system) or HbA1c plus CRP (two systems). This
#' function tabulates the combinations of at-risk biomarkers, their
#' co-occurrence and the redundancy of the panel.
#'
#' @param result An `al_result` from a count protocol (`quantile_count` or
#'   `cutoff_count`).
#' @param top Number of most frequent patterns shown in the plots.
#' @return An `al_patterns` object with `patterns` (each combination of at-risk
#'   biomarkers among complete panels, its load, number of systems, count and
#'   share), `pattern_groups` (shares within groups, when a group column was
#'   scored), `cooccurrence` (pairwise joint risk: observed and
#'   independence-expected counts, their ratio and the phi coefficient),
#'   `redundancy` (eigenvalues of the phi matrix and the effective number of
#'   independent indicators of Li and Ji, 2005) and editable plots:
#'   `frequency`, `matrix`, `patterns` (both combined, requires patchwork) and
#'   `cooccurrence`. Medication-defined risk counts as risk.
#' @references Li J, Ji L (2005). Adjusting multiple testing in multilocus
#'   analyses using the eigenvalues of a correlation matrix. Heredity 95, 221-227.
#' @examples
#' source(system.file("examples", "panel.R", package = "allostaticR"))
#' patterns <- risk_patterns(allostatic_load(data, spec, group = "group"))
#' patterns$patterns
#' patterns$redundancy$meff
#' @family risk patterns
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
risk_patterns <- function(result, top = 15L) {
  .assert_result(result); .check_integer(top, "top", 1L)
  spec <- result$reference$spec
  if (!.is_count_method(spec)) stop("Risk patterns require binary risk flags (quantile_count or cutoff_count).", call. = FALSE)
  if ("group" %in% names(result$scores) && all(is.na(result$scores$group)))
    stop("Group column has no non-missing values.", call. = FALSE)
  info <- .biomarker_info(spec)
  mat <- .risk_matrix(result)
  rows <- which(rowSums(is.na(mat)) == 0L)
  if (length(rows) < 2L) stop("At least two complete panels are required.", call. = FALSE)
  cp <- mat[rows, , drop = FALSE]
  key <- apply(cp, 1, paste0, collapse = "")
  pattern_label <- function(r) if (any(r == 1L)) paste(info$name[r == 1L], collapse = "+") else "none"
  raw_pattern <- vapply(seq_len(nrow(cp)), function(i) pattern_label(cp[i, ]), character(1))
  counts <- table(key)
  first <- match(names(counts), key)
  patterns <- data.frame(
    key = names(counts),
    pattern = make.unique(raw_pattern[first], sep = " #"),
    n = as.integer(counts),
    stringsAsFactors = FALSE
  )
  idx_rows <- match(patterns$key, key)
  patterns$load <- vapply(idx_rows, function(i) sum(cp[i, ] == 1L), integer(1))
  patterns$systems <- vapply(idx_rows, function(i) paste(unique(info$system[cp[i, ] == 1L]), collapse = "+"), character(1))
  patterns$n_systems <- vapply(idx_rows, function(i) length(unique(info$system[cp[i, ] == 1L])), integer(1))
  patterns <- patterns[order(-patterns$n, patterns$load, patterns$pattern), ]
  patterns$share <- patterns$n / length(rows)
  patterns$rank <- seq_len(nrow(patterns))
  pattern_keys <- patterns$key
  pattern_rows <- match(pattern_keys, key)
  patterns <- patterns[c("rank", "pattern", "load", "n_systems", "systems", "n", "share")]
  rownames(patterns) <- NULL
  pattern_groups <- NULL
  if ("group" %in% names(result$scores)) {
    g <- as.character(result$scores$group[rows]); keep <- !is.na(g)
    pattern_groups <- as.data.frame(table(key = factor(key[keep], pattern_keys), group = g[keep]),
      stringsAsFactors = FALSE)
    pattern_groups$rank <- match(pattern_groups$key, pattern_keys)
    pattern_groups$pattern <- patterns$pattern[pattern_groups$rank]
    names(pattern_groups)[3] <- "n"
    pattern_groups <- pattern_groups[c("rank", "pattern", "group", "n")]
    pattern_groups$share <- pattern_groups$n / stats::ave(pattern_groups$n, pattern_groups$group, FUN = sum)
  }
  pairs <- if (ncol(cp) > 1L) utils::combn(ncol(cp), 2) else matrix(integer(), 2, 0)
  cooccurrence <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1, k]; j <- pairs[2, k]; x <- cp[, i]; y <- cp[, j]
    expected <- length(rows) * mean(x) * mean(y)
    data.frame(biomarker_a = info$name[i], biomarker_b = info$name[j],
      same_system = info$system[i] == info$system[j], n_both = sum(x == 1L & y == 1L),
      expected_both = expected, ratio = if (expected > 0) sum(x == 1L & y == 1L) / expected else NA_real_,
      phi = if (stats::sd(x) > 0 && stats::sd(y) > 0) stats::cor(x, y) else NA_real_)
  }))
  variable <- apply(cp, 2, stats::sd) > 0
  redundancy <- list(n_biomarkers = ncol(cp), n_informative = sum(variable), meff = NA_real_, eigenvalues = numeric())
  if (sum(variable) >= 1L) {
    phi <- if (sum(variable) > 1L) stats::cor(cp[, variable, drop = FALSE]) else matrix(1)
    redundancy$eigenvalues <- eigen(phi, symmetric = TRUE, only.values = TRUE)$values
    redundancy$meff <- .meff(redundancy$eigenvalues)
  }
  out <- structure(list(patterns = patterns, pattern_groups = pattern_groups, cooccurrence = cooccurrence,
    redundancy = redundancy, n_complete = length(rows), n_total = nrow(mat), plots = list()), class = "al_patterns")
  out$plots <- .pattern_plots(out, info, top, spec, pattern_rows = pattern_rows, cp = cp)
  if (!is.null(out$pattern_groups)) out$pattern_groups$rank <- NULL
  out
}

.pattern_plots <- function(x, info, top, spec, pattern_rows, cp) {
  shown <- utils::head(x$patterns, top)
  shown_rows <- pattern_rows[shown$rank]
  colors <- .system_colors(.protocol_systems(spec))
  if (!is.null(x$pattern_groups)) {
    bars <- x$pattern_groups[x$pattern_groups$rank %in% shown$rank, ]
    bars$x <- match(bars$rank, shown$rank)
    groups <- sort(unique(bars$group))
    fills <- if (length(groups) <= length(.al_palette$groups)) .al_palette$groups[seq_along(groups)]
      else .palette_for(length(groups), .al_palette$categorical)
    frequency <- ggplot2::ggplot(bars, ggplot2::aes(x = .data$x, y = .data$share, fill = .data$group)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.74) +
      ggplot2::scale_fill_manual(values = stats::setNames(fills, groups))
  } else {
    bars <- shown; bars$x <- seq_len(nrow(shown))
    frequency <- ggplot2::ggplot(bars, ggplot2::aes(x = .data$x, y = .data$share)) +
      ggplot2::geom_col(width = 0.74, fill = .al_palette$secondary)
  }
  frequency_args <- list(title = "Risk patterns behind the load score", x = NULL,
    subtitle = sprintf("Most frequent combinations of at-risk biomarkers; top %d of %d cover %s of %d complete panels",
      nrow(shown), nrow(x$patterns), .percent(sum(shown$n) / x$n_complete), x$n_complete),
    y = if (!is.null(x$pattern_groups)) "Share of complete panels in group" else "Share of complete panels")
  if (!is.null(x$pattern_groups)) frequency_args$fill <- NULL
  frequency <- do.call(.finish_plot, c(list(frequency +
    ggplot2::scale_x_continuous(breaks = NULL, expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(labels = .percent, expand = ggplot2::expansion(mult = c(0, 0.05)))),
    frequency_args))
  membership <- cp[shown_rows, , drop = FALSE] == 1L
  grid <- expand.grid(x = seq_len(nrow(shown)), biomarker = info$name, stringsAsFactors = FALSE)
  grid$member <- membership[cbind(grid$x, match(grid$biomarker, info$name))]
  grid$y <- factor(info$label[match(grid$biomarker, info$name)], levels = rev(info$label))
  grid$system <- info$system[match(grid$biomarker, info$name)]
  links <- do.call(rbind, lapply(seq_len(nrow(shown)), function(k) {
    d <- grid[grid$x == k & grid$member, ]
    if (nrow(d) < 2L) return(NULL)
    data.frame(x = k, from = min(as.integer(d$y)), to = max(as.integer(d$y)))
  }))
  loads <- data.frame(x = seq_len(nrow(shown)), label = paste("AL", shown$load))
  matrix_plot <- .finish_plot(
    ggplot2::ggplot(grid, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(data = grid[!grid$member, ], colour = .al_palette$grid, size = 3) +
      (if (!is.null(links)) ggplot2::geom_segment(data = links, ggplot2::aes(x = .data$x, xend = .data$x,
        y = .data$from, yend = .data$to), inherit.aes = FALSE, colour = .al_palette$secondary, linewidth = 0.6) else NULL) +
      ggplot2::geom_point(data = grid[grid$member, ], ggplot2::aes(colour = .data$system), size = 3.4) +
      ggplot2::geom_text(data = loads, ggplot2::aes(x = .data$x, y = length(info$name) + 0.8, label = .data$label),
        size = 2.8, colour = .al_palette$secondary, inherit.aes = FALSE) +
      ggplot2::scale_colour_manual(values = colors, breaks = names(colors)) +
      ggplot2::scale_x_continuous(breaks = NULL, expand = ggplot2::expansion(add = 0.6)) +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, 1.2))),
    x = NULL, y = NULL, colour = "At risk in pattern", width = 100,
    caption = sprintf("Complete panels only (%d of %d records). Effective number of independent indicators: %.1f of %d.",
      x$n_complete, x$n_total, x$redundancy$meff, x$redundancy$n_biomarkers)) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  plots <- list(frequency = frequency, matrix = matrix_plot)
  if (requireNamespace("patchwork", quietly = TRUE))
    plots$patterns <- patchwork::wrap_plots(frequency + ggplot2::theme(legend.position = "top"), matrix_plot,
      ncol = 1, heights = c(1.15, 1))
  if (!is.null(x$cooccurrence) && nrow(x$cooccurrence)) {
    co <- x$cooccurrence
    co$a <- factor(info$label[match(co$biomarker_a, info$name)], levels = info$label)
    co$b <- factor(info$label[match(co$biomarker_b, info$name)], levels = rev(info$label))
    plots$cooccurrence <- ggplot2::ggplot(co, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$phi)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(.data$phi), "NA", sprintf("%.2f", .data$phi))),
        size = 3, colour = .al_palette$primary) +
      ggplot2::scale_fill_gradient2(low = .al_palette$diverging[["low"]], mid = .al_palette$diverging[["mid"]],
        high = .al_palette$diverging[["high"]], midpoint = 0, limits = c(-1, 1), na.value = .al_palette$unavailable) +
      theme_allostatic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
      ggplot2::labs(title = "Co-occurrence of risk flags", x = NULL, y = NULL, fill = "Phi",
        subtitle = "Phi coefficient between binary risk flags among complete panels",
        caption = .wrap_caption("Positive values: biomarkers at risk together more often than expected under independence.", 70))
  }
  plots
}

#' @export
print.al_patterns <- function(x, ...) {
  cat("Risk patterns among", x$n_complete, "complete panels of", x$n_total, "records\n")
  cat(sprintf("Distinct patterns: %d | effective independent indicators: %.2f of %d\n\n",
    nrow(x$patterns), x$redundancy$meff, x$redundancy$n_biomarkers))
  print(utils::head(x$patterns, 10L), row.names = FALSE, digits = 3, ...)
  invisible(x)
}

#' Map individual risk contributions
#'
#' An oncoprint-style map: one column per participant (sorted by load, then by
#' risk pattern), one row per biomarker grouped by system. Count protocols show
#' at-risk flags coloured by system; z-score protocols show the risk-oriented
#' z-score on a diverging scale. Missing measurements are blank cells. A top
#' panel shows each participant's score and a side panel the prevalence (or
#' mean z) by group.
#'
#' @param result An `al_result`.
#' @param max_participants Maximum number of participant columns. Larger
#'   cohorts are sampled reproducibly within groups and the sample is declared
#'   in the caption; patterns in the map are then illustrative.
#' @param seed Local random seed for that sample.
#' @return An editable patchwork of ggplot objects (requires patchwork), or the
#'   central tile map alone when patchwork is unavailable. Participant
#'   identifiers are not displayed.
#' @family risk patterns
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
plot_risk_map <- function(result, max_participants = 2000L, seed = 1L) {
  .assert_result(result); .check_integer(max_participants, "max_participants", 10L)
  spec <- result$reference$spec; info <- .biomarker_info(spec); count <- .is_count_method(spec)
  scores <- result$scores
  mat <- .contribution_matrix(result)
  group <- if ("group" %in% names(scores)) {
    g <- as.character(scores$group); g[is.na(g)] <- "Missing group"; g
  } else rep("All participants", nrow(scores))
  shown <- seq_len(nrow(scores)); note <- ""
  if (length(shown) > max_participants) {
    shown <- .with_seed(seed, function() unlist(lapply(split(shown, group), function(r)
      r[sort(sample.int(length(r), max(1L, round(max_participants * length(r) / nrow(scores)))))]), use.names = FALSE))
    note <- sprintf(" Showing a random sample of %d of %d participants (seed %d).", length(shown), nrow(scores), seed)
  }
  key <- apply(mat, 1, function(r) { r[is.na(r)] <- -1; sum((if (count) r else r > 0) * 3^(rev(seq_along(r)) - 1)) })
  ord <- shown[order(group[shown], is.na(scores$score[shown]), -ifelse(is.na(scores$score[shown]), 0, scores$score[shown]), -key[shown])]
  position <- integer(nrow(scores))
  for (g in unique(group[ord])) { r <- ord[group[ord] == g]; position[r] <- seq_along(r) }
  tiles <- data.frame(row = rep(shown, times = ncol(mat)), biomarker = rep(info$name, each = length(shown)),
    value = as.vector(mat[shown, , drop = FALSE]))
  tiles <- tiles[!is.na(tiles$value), ]
  tiles$position <- position[tiles$row]; tiles$group <- group[tiles$row]
  tiles$y <- factor(info$label[match(tiles$biomarker, info$name)], levels = rev(info$label))
  tiles$system <- info$system[match(tiles$biomarker, info$name)]
  colors <- .system_colors(.protocol_systems(spec))
  facet <- ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free_x",
    labeller = ggplot2::label_wrap_gen(width = 18))
  if (count) {
    tiles$state <- ifelse(tiles$value == 1, tiles$system, "not at risk")
    fill_scale <- ggplot2::scale_fill_manual(values = c(colors, `not at risk` = .al_palette$neutral),
      breaks = c(names(colors), "not at risk"), labels = c(paste("At risk:", names(colors)), "Not at risk"))
    fill_aes <- ggplot2::aes(fill = .data$state)
  } else {
    tiles$z <- pmin(pmax(tiles$value, -3), 3)
    fill_scale <- ggplot2::scale_fill_gradient2(low = .al_palette$diverging[["low"]], mid = .al_palette$diverging[["mid"]],
      high = .al_palette$diverging[["high"]], midpoint = 0, limits = c(-3, 3), breaks = -3:3)
    fill_aes <- ggplot2::aes(fill = .data$z)
  }
  scored <- data.frame(position = position[shown], group = group[shown], score = scores$score[shown])
  boundary <- do.call(rbind, lapply(split(scored, scored$group), function(d) {
    n_valid <- sum(!is.na(d$score))
    if (n_valid == nrow(d)) return(NULL)
    data.frame(group = d$group[1], x = n_valid + 0.5, mid = (n_valid + nrow(d)) / 2 + 0.5,
      show_label = (nrow(d) - n_valid) >= max(5L, round(0.015 * length(shown))))
  }))
  map <- .finish_plot(
    ggplot2::ggplot(tiles, ggplot2::aes(x = .data$position, y = .data$y)) +
      ggplot2::geom_tile(fill_aes, width = 0.8, height = 0.82) +
      (if (!is.null(boundary)) ggplot2::geom_vline(data = boundary, ggplot2::aes(xintercept = .data$x),
        colour = .al_palette$muted, linewidth = 0.3) else NULL) +
      facet + fill_scale + ggplot2::scale_x_continuous(expand = c(0, 0)),
    x = "Participants sorted by score and risk pattern (blank = missing)",
    y = NULL, fill = if (count) NULL else "Risk-oriented z") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(), panel.spacing.x = ggplot2::unit(7, "pt"))
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("Install patchwork for the score and prevalence panels; returning the tile map only.")
    return(.finish_plot(map, title = "Risk map", width = 100,
      caption = paste0("Reference: ", result$reference$label, ".", note)))
  }
  valid <- scored[!is.na(scored$score), ]
  if (!nrow(valid)) stop("No valid scores to plot.", call. = FALSE)
  n_groups <- length(unique(group))
  top <- .finish_plot(
    ggplot2::ggplot(valid, ggplot2::aes(x = .data$position, y = .data$score)) +
      ggplot2::geom_blank(data = scored, ggplot2::aes(x = .data$position, y = 0)) +
      ggplot2::geom_col(width = 0.8, fill = .al_palette$secondary) +
      (if (!is.null(boundary)) list(ggplot2::geom_vline(data = boundary, ggplot2::aes(xintercept = .data$x),
        colour = .al_palette$muted, linewidth = 0.3),
      if (any(boundary$show_label)) ggplot2::geom_text(
        data = boundary[boundary$show_label, , drop = FALSE],
        ggplot2::aes(x = .data$x, y = max(c(valid$score, 1)) * 0.5, label = "not scored"),
        angle = 90, hjust = 0.5, vjust = -0.45, size = 2.5, colour = .al_palette$muted) else NULL) else NULL) +
      facet + ggplot2::scale_x_continuous(expand = c(0, 0)),
    x = NULL, y = "Score") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
      strip.clip = "off",
      strip.text = if (n_groups >= 3L) {
        ggplot2::element_text(face = "bold", angle = 20, hjust = 0, vjust = 0, size = 8.2)
      } else {
        ggplot2::element_text(face = "bold", hjust = 0.5, size = 8.5)
      },
      panel.spacing.x = ggplot2::unit(9, "pt"))
  side_data <- do.call(rbind, lapply(info$name, function(b) do.call(rbind, lapply(unique(group), function(g) {
    x <- mat[group == g, match(b, info$name)]
    # Inf contributions are reported as NA.
    data.frame(biomarker = b, group = g, value = .mean_or_na(x))
  }))))
  side_data$y <- factor(info$label[match(side_data$biomarker, info$name)], levels = rev(info$label))
  side_data$system <- info$system[match(side_data$biomarker, info$name)]
  groups <- unique(group)
  shapes <- stats::setNames(c(21, 16, 17, 15, 18, 25)[seq_len(min(6L, length(groups)))], utils::head(groups, 6L))
  side_data <- side_data[side_data$group %in% names(shapes) & is.finite(side_data$value), ]
  if (!nrow(side_data)) stop("No valid scores to plot.", call. = FALSE)
  spans <- do.call(rbind, lapply(split(side_data, side_data$biomarker), function(d) data.frame(y = d$y[1],
    system = d$system[1], from = min(d$value, na.rm = TRUE), to = max(d$value, na.rm = TRUE))))
  side <- .finish_plot(
    ggplot2::ggplot(side_data, ggplot2::aes(x = .data$value, y = .data$y)) +
      (if (!count) ggplot2::geom_vline(xintercept = 0, colour = .al_palette$grid) else NULL) +
      ggplot2::geom_segment(data = spans, ggplot2::aes(x = .data$from, xend = .data$to, y = .data$y, yend = .data$y,
        colour = .data$system), linewidth = 1.1) +
      ggplot2::geom_point(ggplot2::aes(shape = .data$group, colour = .data$system), size = 2.4, stroke = 0.9,
        fill = "white", na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = colors, guide = "none") + ggplot2::scale_shape_manual(values = shapes) +
      (if (count) ggplot2::scale_x_continuous(labels = .percent, limits = c(0, max(0.5, side_data$value, na.rm = TRUE)),
        expand = ggplot2::expansion(mult = c(0.04, 0.08)))
        else ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.06)))),
    x = if (count) "At risk (observed)" else "Mean z", y = NULL, shape = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  patchwork::wrap_plots(A = top, B = map, C = side, design = "AAAAA##\nBBBBBCC",
    heights = c(1, 2.2), guides = "collect") +
    patchwork::plot_annotation(title = "Risk map: which biomarkers drive each participant's load",
      subtitle = if (count) "Each column is a participant; the same score can come from different systems"
        else "Each column is a participant; cells show risk-oriented z-scores (capped at \u00b13)",
      caption = .wrap_caption(paste0("Reference: ", result$reference$label, " | ", result$audit$method, " | ",
        result$audit$scale, ".", note), 110),
      theme = theme_allostatic() + ggplot2::theme(legend.position = "bottom", legend.box = "vertical"))
}
