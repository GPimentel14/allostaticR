#' Compare scores for matched participants
#' @param results Named list of at least two al_result objects. Identifier sets
#'   must be identical; order may differ. Cross-scale raw differences are omitted.
#' @return Pairwise comparisons including matched coverage, Spearman correlation
#'   and raw differences only for the same scale and panel size.
#' @family sensitivity analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
compare_scores <- function(results) {
  .named_list(results, 2L); invisible(lapply(results, .assert_result))
  ids <- results[[1]]$scores$id
  for (r in results) if (length(r$scores$id) != length(ids) || !setequal(r$scores$id, ids))
    stop("All results must have the same participant identifier set.", call. = FALSE)
  pairs <- utils::combn(seq_along(results), 2)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    a <- results[[pairs[1, j]]]; b <- results[[pairs[2, j]]]
    x <- a$scores$score[match(ids, a$scores$id)]; y <- b$scores$score[match(ids, b$scores$id)]
    ok <- is.finite(x) & is.finite(y)
    same_scale <- identical(a$audit$scale, b$audit$scale) &&
      (a$audit$scale != "count" || a$scores$n_total[1] == b$scores$n_total[1])
    data.frame(result_a = names(results)[pairs[1, j]], result_b = names(results)[pairs[2, j]],
      n_matched = sum(ok), n_valid_a = sum(is.finite(x)), n_valid_b = sum(is.finite(y)),
      spearman = .safe_spearman(x, y), same_scale = same_scale,
      mean_difference_b_minus_a = if (same_scale) .mean_or_na(y[ok] - x[ok]) else NA_real_,
      mean_absolute_difference = if (same_scale) .mean_or_na(abs(y[ok] - x[ok])) else NA_real_)
  }))
}

#' Evaluate sensitivity to alternative protocols
#' @param data Study data frame.
#' @param reference_data Independent calibration data frame.
#' @param specs Named list of at least two define_protocol objects. The first is baseline.
#' @param id Optional participant identifier column.
#' @return An `al_sensitivity` object containing fitted results, matched comparisons,
#'   normalized changes and editable plots. Counts are divided by panel size;
#'   z-score scales are not risk fractions and have no normalized change.
#'   `plots$changes` orders participants by baseline risk fraction;
#'   `plots$reclassification` cross-tabulates count scenarios against a count baseline.
#' @family sensitivity analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
sensitivity_analysis <- function(data, reference_data, specs, id = NULL) {
  .named_list(specs, 2L); invisible(lapply(specs, validate_protocol))
  results <- lapply(specs, function(s) score_load(data, fit_reference(reference_data, s), id = id))
  base <- .risk_fraction(results[[1]])
  changes <- do.call(rbind, lapply(seq_along(results), function(i) {
    x <- .risk_fraction(results[[i]])
    data.frame(row = results[[i]]$scores$row, scenario = names(results)[i], normalized_score = x, change = x - base)
  }))
  plots <- list()
  if (any(is.finite(changes$change))) {
    order_rank <- rank(ifelse(is.finite(base), base, Inf), ties.method = "first")
    changes_plot <- changes
    changes_plot$position <- order_rank[changes_plot$row]
    changes_plot$scenario <- factor(changes_plot$scenario, levels = rev(names(results)))
    plots$changes <- ggplot2::ggplot(changes_plot,
        ggplot2::aes(x = .data$position, y = .data$scenario, fill = .data$change)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = .al_palette$diverging[["low"]], mid = .al_palette$diverging[["mid"]],
        high = .al_palette$diverging[["high"]], midpoint = 0, na.value = .al_palette$unavailable,
        limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      theme_allostatic() + ggplot2::labs(title = "Protocol sensitivity",
        x = "Participants ordered by baseline risk fraction (unscored last)",
        y = "Protocol scenario", fill = "Change in\nrisk fraction",
        caption = .wrap_caption(paste0("Each scenario minus '", names(results)[1],
          "'. Light gray: no change; dark gray: unavailable paired score (missing panel or z-score scale).")))
  }
  if (results[[1]]$audit$scale == "count") {
    tabs <- do.call(rbind, lapply(names(results)[-1], function(nm) {
      r <- results[[nm]]
      if (r$audit$scale != "count") return(NULL)
      x <- results[[1]]$scores$score; y <- r$scores$score; ok <- is.finite(x) & is.finite(y)
      if (!any(ok)) return(NULL)
      tab <- as.data.frame(table(baseline = factor(x[ok], 0:results[[1]]$scores$n_total[1]),
        scenario_score = factor(y[ok], 0:r$scores$n_total[1])), stringsAsFactors = FALSE)
      tab$scenario <- nm
      tab
    }))
    if (!is.null(tabs) && nrow(tabs)) {
      tabs$baseline <- as.integer(tabs$baseline); tabs$scenario_score <- as.integer(tabs$scenario_score)
      tabs$diagonal <- tabs$baseline == tabs$scenario_score
      plots$reclassification <- ggplot2::ggplot(tabs[tabs$Freq > 0, ],
          ggplot2::aes(x = .data$baseline, y = .data$scenario_score)) +
        ggplot2::geom_tile(ggplot2::aes(fill = .data$diagonal), width = 0.92, height = 0.92) +
        ggplot2::geom_text(ggplot2::aes(label = .data$Freq), size = 3.2, colour = .al_palette$primary) +
        ggplot2::scale_fill_manual(values = c(`TRUE` = .al_palette$light_blue, `FALSE` = .al_palette$light_red),
          labels = c(`TRUE` = "Same score", `FALSE` = "Reclassified")) +
        ggplot2::facet_wrap(~scenario) + ggplot2::coord_equal() + theme_allostatic() +
        ggplot2::labs(title = "Reclassification between protocols", x = paste("Baseline count:", names(results)[1]),
          y = "Scenario count", fill = NULL, caption = "Numbers: participants with valid scores under both protocols.")
    }
  }
  structure(list(results = results, comparisons = compare_scores(results), changes = changes,
    plots = plots, baseline = names(specs)[1]), class = "al_sensitivity")
}

#' Measure influence by leaving out biomarkers or systems
#' @param data Study data frame.
#' @param reference Fixed fitted reference.
#' @param by Omit one biomarker or one system at a time.
#' @param id Optional participant identifier column.
#' @return An `al_score_comparison` object with paired comparisons and an editable
#'   influence plot. Count scores are compared as risk fractions; other scales
#'   in their own units. Thresholds are never refitted.
#' @family sensitivity analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
panel_influence <- function(data, reference, by = c("biomarker", "system"), id = NULL) {
  .assert_reference(reference); by <- match.arg(by)
  key <- .biomarker_field(reference$spec, if (by == "biomarker") "name" else "system")
  units <- unique(key)
  if (length(units) < 2L) stop("At least two biomarkers/systems required.", call. = FALSE)
  base <- score_load(data, reference, id = id)
  x <- .comparable_score(base)
  results <- lapply(units, function(u) score_load(data, .subset_reference(reference, key != u), id = id))
  names(results) <- units
  summary <- do.call(rbind, lapply(seq_along(results), function(i) {
    y <- .comparable_score(results[[i]]); ok <- is.finite(x) & is.finite(y)
    data.frame(omitted = units[i], n_matched = sum(ok), n_valid_without = sum(is.finite(y)),
      mean_change = .mean_or_na(y[ok] - x[ok]), mean_absolute_change = .mean_or_na(abs(y[ok] - x[ok])),
      spearman = .safe_spearman(x, y))
  }))
  changes <- do.call(rbind, lapply(seq_along(results), function(i)
    data.frame(row = base$scores$row, id = base$scores$id, omitted = units[i],
      baseline_fraction = x, reduced_fraction = .comparable_score(results[[i]]))))
  p <- summary[is.finite(summary$mean_absolute_change), ]
  unit_label <- if (base$audit$scale == "count") "risk fraction" else .scale_label(base$audit$scale)
  plots <- list()
  if (nrow(p)) plots$influence <- ggplot2::ggplot(p,
      ggplot2::aes(x = .data$mean_absolute_change, y = .data$omitted)) +
    ggplot2::geom_point(color = .al_palette$groups[["ink"]], size = 3) + theme_allostatic() +
    ggplot2::labs(title = "Panel influence", x = paste("Mean absolute change in", tolower(unit_label)),
      y = paste("Omitted", by), caption = "Fixed reference thresholds; paired valid rows only. Denominator changes with the panel.")
  structure(list(baseline = base, results = results, summary = summary, changes = changes,
    plots = plots, by = by), class = "al_influence")
}
