#' Harmonize biomarker distributions across cohorts
#'
#' Maps each biomarker of every cohort onto the distribution of a target
#' cohort before pooling, for multi-cohort and individual-participant
#' meta-analyses whose laboratories or protocols differ. `"quantile"` equating
#' replaces each value by the target quantile at its within-cohort mid-rank;
#' `"mean_sd"` matches the mean and SD on the transformed scale. Harmonization
#' removes all distributional differences between cohorts, including real
#' population differences, so use it only when cohorts are assumed
#' exchangeable (for example random splits of one population measured in
#' different laboratories), and report results with and without it.
#'
#' @param data Pooled data frame with a cohort column.
#' @param spec Protocol defining the biomarkers (and transforms).
#' @param cohort Cohort column.
#' @param target Target cohort level; defaults to the largest cohort.
#' @param method `"quantile"` or `"mean_sd"`.
#' @param strata Optional columns within which harmonization is performed
#'   separately (for example sex).
#' @return An `al_harmonization` with `data` (harmonized copy; original
#'   values kept as `<biomarker>_original`), `diagnostics` (per biomarker and
#'   cohort: medians and Kolmogorov-Smirnov distance to the target before and
#'   after) and an editable plot of distributions before and after.
#' @family multi-cohort analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
harmonize_cohorts <- function(data, spec, cohort, target = NULL, method = c("quantile", "mean_sd"), strata = character()) {
  validate_protocol(spec); method <- match.arg(method); .scalar_string(cohort, "cohort")
  if (!is.data.frame(data) || !all(c(cohort, strata) %in% names(data))) stop("Missing cohort or strata columns.", call. = FALSE)
  if (anyNA(data[[cohort]]) || any(vapply(strata, function(s) anyNA(data[[s]]), logical(1)))) stop("Cohort and strata columns must be complete.", call. = FALSE)
  validate_data(data, define_protocol(spec$biomarkers))
  cohorts <- as.character(data[[cohort]])
  if (length(unique(cohorts)) < 2L) stop("At least two cohorts are required.", call. = FALSE)
  if (is.null(target)) target <- names(sort(table(cohorts), decreasing = TRUE))[1]
  if (!target %in% cohorts) stop("target must be an existing cohort.", call. = FALSE)
  key <- .stratum_key(data, strata)
  out <- data; diagnostics <- list(); long <- list()
  ks <- function(a, b) if (length(a) > 1L && length(b) > 1L) unname(suppressWarnings(stats::ks.test(a, b)$statistic)) else NA_real_
  for (b in spec$biomarkers) {
    x <- data[[b$name]]; new <- x
    for (k in unique(key)) {
      ref <- x[cohorts == target & key == k]; ref <- ref[!is.na(ref)]
      if (length(ref) < 2L) stop("Too few target values for ", b$name, " in stratum ", k, ".", call. = FALSE)
      for (coh in setdiff(unique(cohorts[key == k]), target)) {
        rows <- which(cohorts == coh & key == k & !is.na(x))
        if (length(rows) < 2L) next
        if (method == "quantile") {
          p <- (rank(x[rows]) - 0.5) / length(rows)
          new[rows] <- unname(stats::quantile(ref, p, type = 7))
        } else {
          tr <- .transform(ref, b$transform); tx <- .transform(x[rows], b$transform)
          sd_tr <- stats::sd(tr); sd_tx <- stats::sd(tx)
          if (!is.finite(sd_tr) || sd_tr <= 0) stop(sprintf("Target cohort '%s' has no variance for biomarker '%s'.", target, b$name), call. = FALSE)
          if (!is.finite(sd_tx) || sd_tx <= 0) stop(sprintf("Cohort '%s' has no variance for biomarker '%s'.", coh, b$name), call. = FALSE)
          new[rows] <- .untransform((tx - mean(tx)) / sd_tx * sd_tr + mean(tr), b$transform)
        }
      }
    }
    out[[paste0(b$name, "_original")]] <- x; out[[b$name]] <- new
    for (coh in unique(cohorts)) {
      r <- x[cohorts == target]; before <- x[cohorts == coh]; after <- new[cohorts == coh]
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(biomarker = b$name, cohort = coh, is_target = coh == target,
        n = sum(!is.na(before)), median_before = stats::median(before, na.rm = TRUE), median_after = stats::median(after, na.rm = TRUE),
        median_target = stats::median(r, na.rm = TRUE), ks_before = if (coh == target) 0 else ks(stats::na.omit(before), stats::na.omit(r)),
        ks_after = if (coh == target) 0 else ks(stats::na.omit(after), stats::na.omit(r)))
      b_vals <- as.numeric(stats::na.omit(before)); a_vals <- as.numeric(stats::na.omit(after))
      if (length(b_vals)) long[[length(long) + 1L]] <- data.frame(biomarker = b$label, cohort = coh, stage = "Before", value = b_vals)
      if (length(a_vals)) long[[length(long) + 1L]] <- data.frame(biomarker = b$label, cohort = coh, stage = "After", value = a_vals)
    }
  }
  attr(out, "al_harmonization") <- list(method = method, target = target, cohort = cohort, strata = strata)
  long <- if (length(long)) do.call(rbind, long) else data.frame(biomarker = character(), cohort = character(), stage = character(), value = numeric())
  long$stage <- factor(long$stage, c("Before", "After"))
  cohort_levels <- c(target, setdiff(unique(cohorts), target))
  colors <- stats::setNames(c(.al_palette$secondary, .palette_for(length(cohort_levels) - 1L, .al_palette$categorical)), cohort_levels)
  plot <- .finish_plot(
    ggplot2::ggplot(long, ggplot2::aes(x = .data$value, colour = .data$cohort)) +
      ggplot2::geom_density(linewidth = 0.7, na.rm = TRUE) +
      ggplot2::facet_grid(stage ~ biomarker, scales = "free") + ggplot2::scale_colour_manual(values = colors),
    title = "Biomarker distributions before and after harmonization", x = "Value (declared unit)",
    y = "Density", colour = "Cohort", width = 100,
    caption = sprintf("Method: %s; target cohort: %s%s. Harmonization also removes real population differences.",
      method, target, if (length(strata)) paste0("; within ", paste(strata, collapse = ", ")) else ""))
  d <- do.call(rbind, diagnostics); rownames(d) <- NULL
  structure(list(data = out, diagnostics = d, method = method, target = target, plots = list(distributions = plot)),
    class = "al_harmonization")
}

#' @export
print.al_harmonization <- function(x, ...) {
  cat(sprintf("Cohort harmonization (%s) to target cohort '%s'\n\n", x$method, x$target))
  print(x$diagnostics, row.names = FALSE, digits = 3, ...)
  invisible(x)
}
