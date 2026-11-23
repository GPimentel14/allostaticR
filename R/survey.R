#' Survey-weighted summaries of allostatic load
#'
#' Population estimates from complex surveys (for example NHANES): weighted
#' mean score and, for count protocols, weighted prevalence of each biomarker
#' at risk, overall and by group, with percentile intervals from a cluster
#' bootstrap that resamples primary sampling units within design strata.
#'
#' @param result An `al_result` scored from `data`.
#' @param data The scored data frame (same rows and order as the result),
#'   containing the design columns.
#' @param weights Column of positive sampling weights.
#' @param design_strata,cluster Optional design strata and primary sampling
#'   unit columns. Without `cluster`, rows are resampled within strata.
#' @param by Optional grouping column in `data`.
#' @param B Bootstrap resamples (0 skips intervals).
#' @param conf Confidence level.
#' @inheritParams simulate_biomarkers
#' @return An `al_survey` object with `estimates` (scope, target, group,
#'   estimate, interval, unweighted n and weighted total) and an editable plot.
#'   `scope` is `"total"` for the overall estimate and `"group"` for each level
#'   of `by`, so a group named `"All participants"` remains distinct from the
#'   total. Only participants with valid scores contribute to score means;
#'   prevalences use observed measurements. The bootstrap omits
#'   finite-population corrections and weight recalibration, so intervals
#'   approximate design-based ones.
#' @family survey analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
survey_summary <- function(result, data, weights, design_strata = NULL, cluster = NULL, by = NULL,
                           B = 500L, conf = 0.95, seed = 1L) {
  .assert_result(result); .check_integer(B, "B", 0L); .check_open_unit(conf, "conf")
  if (!is.data.frame(data) || nrow(data) != nrow(result$scores)) stop("data must be the scored data frame.", call. = FALSE)
  .scalar_string(weights, "weights")
  for (column in c(weights, design_strata, cluster, by)) {
    if (!column %in% names(data)) stop("Missing column: ", column, call. = FALSE)
    if (anyNA(data[[column]])) stop("Design and group columns must be complete: ", column, call. = FALSE)
  }
  w <- data[[weights]]
  if (!is.numeric(w) || any(!is.finite(w)) || any(w <= 0)) stop("Weights must be positive and finite.", call. = FALSE)
  spec <- result$reference$spec; count <- .is_count_method(spec)
  group <- if (is.null(by)) rep(NA_character_, nrow(data)) else as.character(data[[by]])
  grp_levels <- if (is.null(by)) character() else sort(unique(group))
  groups_meta <- data.frame(
    gid = seq_len(1L + length(grp_levels)),
    label = c("All participants", grp_levels),
    is_total = c(TRUE, rep(FALSE, length(grp_levels))),
    stringsAsFactors = FALSE
  )
  score <- result$scores$score; risk <- if (count) .risk_matrix(result) else NULL
  names_all <- .biomarker_field(spec, "name")
  estimate <- function(rows) {
    out <- list()
    for (idx in seq_len(nrow(groups_meta))) {
      in_g <- if (groups_meta$is_total[idx]) rows else rows[group[rows] == groups_meta$label[idx]]
      ok <- in_g[!is.na(score[in_g])]
      out[[length(out) + 1L]] <- data.frame(gid = groups_meta$gid[idx], target = "mean score", group = groups_meta$label[idx],
        value = if (length(ok)) stats::weighted.mean(score[ok], w[ok]) else NA_real_, stringsAsFactors = FALSE)
      if (count) for (j in seq_along(names_all)) {
        okj <- in_g[!is.na(risk[in_g, j])]
        out[[length(out) + 1L]] <- data.frame(gid = groups_meta$gid[idx], target = paste("prevalence:", names_all[j]),
          group = groups_meta$label[idx],
          value = if (length(okj)) stats::weighted.mean(risk[okj, j], w[okj]) else NA_real_, stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, out)
  }
  all_rows <- seq_len(nrow(data))
  point <- estimate(all_rows)
  stratum <- if (is.null(design_strata)) rep("all", nrow(data)) else as.character(data[[design_strata]])
  plan <- .cluster_plan(nrow(data), stratum, if (is.null(cluster)) NULL else data[[cluster]])
  draws <- if (B > 0) .with_seed(seed, function() vapply(seq_len(B), function(i) {
    rows <- .cluster_draw(plan, cluster)
    e <- estimate(rows); e$value[match(paste(point$target, point$gid), paste(e$target, e$gid))]
  }, numeric(nrow(point)))) else NULL
  ci <- if (is.null(draws)) matrix(NA_real_, 2, nrow(point)) else apply(matrix(draws, nrow = nrow(point)), 1, .percentile_interval, conf = conf)
  n_info <- t(vapply(seq_len(nrow(point)), function(i) {
    idx <- point$gid[i]
    rows <- if (groups_meta$is_total[idx]) all_rows else all_rows[group == groups_meta$label[idx]]
    rows <- if (point$target[i] == "mean score") rows[!is.na(score[rows])] else rows[!is.na(risk[rows, match(sub("prevalence: ", "", point$target[i]), names_all)])]
    c(length(rows), sum(w[rows]))
  }, numeric(2)))
  estimates <- data.frame(scope = ifelse(groups_meta$is_total[point$gid], "total", "group"),
    target = point$target, group = point$group, estimate = point$value, lower = ci[1, ],
    upper = ci[2, ], n = n_info[, 1], weighted_n = n_info[, 2])
  shown <- estimates[estimates$target != "mean score", ]
  shown$series <- ifelse(shown$scope == "total", "Total", paste0("Group: ", shown$group))
  plot <- if (nrow(shown)) .finish_plot(
    ggplot2::ggplot(shown, ggplot2::aes(x = .data$estimate, y = sub("prevalence: ", "", .data$target),
      colour = .data$series, xmin = .data$lower, xmax = .data$upper)) +
      ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5), na.rm = TRUE) +
      ggplot2::scale_x_continuous(labels = .percent) +
      ggplot2::scale_colour_manual(values = rep(c(.al_palette$secondary, .al_palette$categorical), length.out = length(unique(shown$series)))),
    title = "Survey-weighted prevalence of biomarker risk", x = "Weighted prevalence",
    y = NULL, colour = NULL, width = 100,
    caption = sprintf("Weights: %s%s%s. %g%% cluster-bootstrap intervals (B = %d).",
      weights, if (!is.null(design_strata)) paste(", strata:", design_strata) else "",
      if (!is.null(cluster)) paste(", PSU:", cluster) else "", 100 * conf, B))
  else NULL
  structure(list(estimates = estimates, design = list(weights = weights, design_strata = design_strata, cluster = cluster, B = B),
    plots = list(prevalence = plot)), class = "al_survey")
}

#' @export
print.al_survey <- function(x, ...) {
  cat("Survey-weighted allostatic load summary (weights:", x$design$weights, ")\n\n")
  print(x$estimates, row.names = FALSE, digits = 3, ...)
  invisible(x)
}
