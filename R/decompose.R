# Weights that turn complete-panel contributions into the score.
.score_weights <- function(spec, scale) {
  systems <- .biomarker_field(spec, "system")
  n <- length(systems)
  switch(scale, count = rep(1, n), proportion = , z = rep(1 / n, n),
    system_proportion = , system_z = 1 / (length(unique(systems)) * as.numeric(table(systems)[systems])))
}

.contribution_matrix <- function(result) {
  matrix(result$components$contribution, nrow = nrow(result$scores))
}

#' Decompose a group difference in allostatic load
#'
#' With complete panels, every score in the package is a weighted sum of
#' biomarker contributions (risk flags or risk-oriented z-scores). The
#' difference in mean score between two groups therefore splits exactly into
#' biomarker or system contributions: weight times the difference in mean
#' contribution. For counts, each contribution is a difference in risk
#' prevalence. The identity is verified before results are returned.
#'
#' @param result An `al_result` scored with a group column.
#' @param levels Two group levels; contributions describe `levels[2]` minus
#'   `levels[1]`. Defaults to the two groups present, in sorted order.
#' @param by Contributions per `"biomarker"` or per `"system"`.
#' @param B Bootstrap resamples within each group for percentile intervals;
#'   0 skips intervals.
#' @param conf Confidence level of the intervals.
#' @inheritParams simulate_biomarkers
#' @return An `al_decomposition` with `summary` (group means, difference and
#'   interval), `contributions` (per term: weight, group means of the
#'   contribution, contribution, interval and share of the gap) and an
#'   editable waterfall plot in `plots$waterfall`. Only participants with
#'   complete panels enter the decomposition; partial panels scored under a
#'   proportion policy are counted in `n_excluded`. The decomposition is
#'   descriptive and does not adjust for confounding.
#' @examples
#' source(system.file("examples", "panel.R", package = "allostaticR"))
#' result <- allostatic_load(data, spec, id = "id", group = "group")
#' gap <- decompose_gap(result, B = 200)
#' gap$contributions
#' gap$plots$waterfall
#' @family group comparison
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
decompose_gap <- function(result, levels = NULL, by = c("biomarker", "system"), B = 1000L, conf = 0.95, seed = 1L) {
  .assert_result(result); by <- match.arg(by)
  levels <- .group_levels(result, levels)
  .check_integer(B, "B", min = 0L); .check_open_unit(conf, "conf")
  spec <- result$reference$spec; scores <- result$scores
  names_all <- .biomarker_field(spec, "name")
  systems <- .biomarker_field(spec, "system")
  group <- as.character(scores$group)
  complete <- !is.na(scores$score) & scores$n_observed == scores$n_total
  rows_a <- which(complete & group %in% levels[1]); rows_b <- which(complete & group %in% levels[2])
  if (length(rows_a) < 2L || length(rows_b) < 2L)
    stop("Each group needs at least two participants with complete panels.", call. = FALSE)
  mat <- .contribution_matrix(result); w <- .score_weights(spec, result$audit$scale)
  used <- c(rows_a, rows_b)
  if (any(abs(drop(mat[used, , drop = FALSE] %*% w) - scores$score[used]) > 1e-8))
    stop("Scores are not additive in their contributions.", call. = FALSE)
  term <- if (by == "biomarker") names_all else systems
  terms <- unique(term)
  gap <- function(a, b) {
    d <- w * (colMeans(mat[b, , drop = FALSE]) - colMeans(mat[a, , drop = FALSE]))
    vapply(terms, function(t) sum(d[term == t]), numeric(1))
  }
  estimate <- gap(rows_a, rows_b)
  draws <- if (B > 0) .with_seed(seed, function() matrix(vapply(seq_len(B), function(i)
    gap(rows_a[sample.int(length(rows_a), replace = TRUE)], rows_b[sample.int(length(rows_b), replace = TRUE)]),
    numeric(length(terms))), nrow = length(terms))) else NULL
  interval <- function(x) if (is.null(draws)) c(NA_real_, NA_real_) else .percentile_interval(x, conf)
  part <- function(rows) vapply(terms, function(t) sum((w * colMeans(mat[rows, , drop = FALSE]))[term == t]), numeric(1))
  component_mean <- function(rows) if (by == "biomarker") colMeans(mat[rows, , drop = FALSE]) else part(rows)
  difference <- mean(scores$score[rows_b]) - mean(scores$score[rows_a])
  total <- interval(if (is.null(draws)) NULL else colSums(draws))
  summary <- data.frame(from = levels[1], to = levels[2], n_from = length(rows_a), n_to = length(rows_b),
    n_excluded = sum(!is.na(scores$score) & group %in% levels & !complete),
    mean_from = mean(scores$score[rows_a]), mean_to = mean(scores$score[rows_b]),
    difference = difference, lower = total[1], upper = total[2], conf = conf, B = B,
    scale = result$audit$scale)
  ci <- if (is.null(draws)) matrix(NA_real_, 2L, length(terms)) else apply(draws, 1, .percentile_interval, conf = conf)
  labels <- .biomarker_field(spec, "label")
  contributions <- data.frame(term = terms,
    label = if (by == "biomarker") labels else terms,
    system = if (by == "biomarker") systems else terms,
    weight = if (by == "biomarker") w else NA_real_,
    mean_from = unname(component_mean(rows_a)), mean_to = unname(component_mean(rows_b)),
    contribution = unname(estimate), lower = ci[1, ], upper = ci[2, ],
    share = if (abs(difference) > 0) unname(estimate) / difference else NA_real_)
  rownames(contributions) <- NULL
  structure(list(summary = summary, contributions = contributions, by = by,
    plots = list(waterfall = .waterfall_plot(summary, contributions, spec))), class = "al_decomposition")
}

.waterfall_plot <- function(summary, contributions, spec) {
  k <- nrow(contributions)
  ends <- summary$mean_from + cumsum(contributions$contribution)
  steps <- data.frame(x = seq_len(k + 2L), label = c(summary$from, contributions$label, summary$to),
    kind = c("total", contributions$system, "total"),
    start = c(0, summary$mean_from, ends[-k], 0), end = c(summary$mean_from, ends, summary$mean_to))
  steps$lo <- c(NA, steps$start[2:(k + 1L)] + contributions$lower, NA)
  steps$hi <- c(NA, steps$start[2:(k + 1L)] + contributions$upper, NA)
  steps$text <- ifelse(steps$kind == "total", sprintf("%.2f", steps$end), sprintf("%+.2f", steps$end - steps$start))
  span <- diff(range(c(0, steps$start, steps$end, steps$lo, steps$hi), na.rm = TRUE))
  steps$text_y <- pmax(steps$start, steps$end, steps$hi, na.rm = TRUE) + 0.04 * span
  connect <- data.frame(x = steps$x[seq_len(k + 1L)] + 0.36, xend = steps$x[seq_len(k + 1L) + 1L] - 0.36,
    y = steps$end[seq_len(k + 1L)])
  colors <- .system_colors(.protocol_systems(spec))
  ggplot2::ggplot(steps) +
    ggplot2::geom_hline(yintercept = 0, colour = .al_palette$grid, linewidth = 0.4) +
    ggplot2::geom_segment(data = connect, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      colour = .al_palette$muted, linewidth = 0.3) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$x - 0.36, xmax = .data$x + 0.36,
      ymin = pmin(.data$start, .data$end), ymax = pmax(.data$start, .data$end), fill = .data$kind)) +
    ggplot2::geom_errorbar(ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi), width = 0.14,
      colour = .al_palette$primary, linewidth = 0.4, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$text_y, label = .data$text), size = 3.2,
      colour = .al_palette$primary) +
    ggplot2::scale_fill_manual(values = c(colors, total = .al_palette$secondary), breaks = names(colors)) +
    ggplot2::scale_x_continuous(breaks = steps$x, labels = function(x) vapply(steps$label[x], function(y)
      paste(strwrap(y, 14), collapse = "\n"), character(1)), expand = ggplot2::expansion(add = 0.5)) +
    theme_allostatic() +
    ggplot2::labs(title = sprintf("Decomposition of the difference: %s minus %s", summary$to, summary$from),
      subtitle = if (is.finite(summary$lower)) sprintf("Mean %.2f \u2192 %.2f; difference %+.2f (%g%% CI %.2f to %.2f) = sum of contributions",
          summary$mean_from, summary$mean_to, summary$difference, 100 * summary$conf, summary$lower, summary$upper)
        else sprintf("Mean %.2f \u2192 %.2f; difference %+.2f = sum of contributions", summary$mean_from, summary$mean_to, summary$difference),
      x = NULL, y = paste("Mean score:", tolower(.scale_label(summary$scale))), fill = "System",
      caption = .wrap_caption(sprintf(paste0("Complete panels (n = %d and %d). Error bars: %g%% percentile intervals from %d ",
        "bootstrap resamples within group. Descriptive; not adjusted for confounding."),
        summary$n_from, summary$n_to, 100 * summary$conf, summary$B), 100))
}

#' @export
print.al_decomposition <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Decomposition of %s minus %s (%s scale): %+.3f", s$to, s$from, s$scale, s$difference))
  if (is.finite(s$lower)) cat(sprintf(" [%g%% CI %.3f, %.3f]", 100 * s$conf, s$lower, s$upper))
  cat(sprintf("\nComplete panels: %d and %d; excluded partial panels: %d\n\n", s$n_from, s$n_to, s$n_excluded))
  print(x$contributions, row.names = FALSE, digits = 3, ...)
  invisible(x)
}
