#' Bootstrap uncertainty in reference calibration
#' @param data Fixed evaluation data frame independent of reference_data.
#' @param reference_data Calibration data with independent rows.
#' @param spec Quantile or z-score protocol from define_protocol().
#' @param B Number of bootstrap resamples, at least two.
#' @param conf Confidence level for percentile intervals.
#' @param seed Local random seed, preserving the caller's RNG state.
#' @param cluster Optional column of primary sampling units in `reference_data`;
#'   whole clusters are resampled (with their weights) instead of rows.
#' @param design_strata Optional column of survey design strata; resampling
#'   happens within them (protocol strata are used when omitted).
#' @return An bootstrap_load object with reference-parameter intervals, resampled
#'   cohort means, failures and an editable plot. These quantify calibration
#'   uncertainty conditional on the fixed study sample, not laboratory error,
#'   individual clinical uncertainty or total sampling uncertainty (see
#'   [measurement_uncertainty()] for biological and analytical variation).
#'   Invalid resamples are reported, never silently replaced. Rows are
#'   resampled within protocol strata, or clusters within design strata when
#'   `cluster` is given (a with-replacement cluster bootstrap without
#'   finite-population correction). Parameter intervals are omitted for
#'   age-conditional references.
#' @family calibration uncertainty
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
bootstrap_load <- function(data, reference_data, spec, B = 200L, conf = 0.95, seed = 1L,
                           cluster = NULL, design_strata = NULL) {
  validate_protocol(spec)
  if (spec$method == "cutoff_count") stop("Bootstrap calibration requires a quantile or z-score protocol.", call. = FALSE)
  validate_data(data, spec); validate_data(reference_data, spec)
  .check_integer(B, "B", 2L, text = "B must be an integer >= 2.")
  .check_open_unit(conf, "conf")
  original <- fit_reference(reference_data, spec)
  estimate <- .mean_or_na(score_load(data, original)$scores$score)
  for (column in c(cluster, design_strata)) {
    .scalar_string(column, "design column")
    if (!column %in% names(reference_data) || anyNA(reference_data[[column]])) stop("Design columns must be complete: ", column, call. = FALSE)
  }
  stratum <- if (is.null(design_strata)) .stratum_key(reference_data, spec$strata) else as.character(reference_data[[design_strata]])
  plan <- .cluster_plan(nrow(reference_data), stratum, if (is.null(cluster)) NULL else reference_data[[cluster]])
  draws <- .with_seed(seed, function() lapply(seq_len(B), function(i) {
    rows <- .cluster_draw(plan, cluster)
    sampled <- reference_data[rows, , drop = FALSE]
    tryCatch({
      ref <- fit_reference(sampled, spec)
      result <- score_load(data, ref)
      list(mean = .mean_or_na(result$scores$score), cutoffs = ref$cutoffs, error = NA_character_)
    }, error = function(e) list(mean = NA_real_, cutoffs = NULL, error = conditionMessage(e)))
  }))
  means <- vapply(draws, `[[`, numeric(1), "mean")
  failures <- data.frame(replicate = seq_len(B), error = vapply(draws, `[[`, character(1), "error"))
  failures <- failures[!is.na(failures$error), , drop = FALSE]
  valid <- Filter(function(x) !is.null(x$cutoffs), draws)
  cutoff_intervals <- if (length(spec$age)) data.frame() else do.call(rbind, lapply(seq_len(nrow(original$cutoffs)), function(j) {
    d <- original$cutoffs[j, ]
    parameters <- if (spec$method == "z_mean") c("center", "scale")
      else if (d$direction == "both") c("lower", "upper") else if (d$direction == "low") "lower" else "upper"
    do.call(rbind, lapply(parameters, function(parameter) {
      values <- vapply(valid, function(x) x$cutoffs[[parameter]][j], numeric(1))
      ci <- .percentile_interval(values, conf)
      out <- data.frame(biomarker = d$biomarker, tail = parameter, unit = d$unit, estimate = d[[parameter]],
        lower = ci[1], upper = ci[2], n_valid = length(values))
      if (length(spec$strata)) out$stratum <- d$stratum
      out
    }))
  }))
  ci <- .percentile_interval(means, conf)
  summary <- data.frame(estimate = estimate, lower = ci[1], upper = ci[2], conf = conf,
    n_valid = sum(is.finite(means)), n_failed_reference = nrow(failures),
    n_undefined_mean = sum(is.na(means)) - nrow(failures), B = B)
  plots <- list()
  if (is.finite(estimate) && all(is.finite(ci))) plots$mean_interval <- .finish_plot(
    ggplot2::ggplot(summary,
      ggplot2::aes(x = "Fixed study cohort", y = .data$estimate, ymin = .data$lower, ymax = .data$upper)) +
      ggplot2::geom_pointrange(color = .al_palette$groups[["ink"]]),
    title = "Reference calibration uncertainty", x = "Evaluation sample",
    y = paste("Mean score:", .scale_label(.score_scale(spec))),
    caption = paste0(conf * 100, "% percentile interval; independent reference-row bootstrap. Valid draws: ", sum(is.finite(means)), "/", B, "."))
  structure(list(summary = summary, cutoff_intervals = cutoff_intervals,
    draws = data.frame(replicate = seq_len(B), mean_score = means), failures = failures,
    reference = original, seed = seed, design = list(cluster = cluster, design_strata = design_strata),
    plots = plots), class = "al_bootstrap")
}
