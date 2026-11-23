#' Detect biomarker outliers and their impact on allostatic load
#'
#' Flags univariate outliers per biomarker (Tukey fences or robust z-scores
#' from the median and MAD, on the protocol transform scale, within protocol
#' strata) and multivariate outliers of the whole panel (robust Mahalanobis
#' distance from a minimum covariance determinant fit, compared with a
#' chi-square cut-off). Nothing is removed or changed unless `winsorize = TRUE`,
#' which returns an explicit copy with values clipped to the fences. With a
#' fitted reference, the impact of winsorizing on each score is reported.
#' Outliers mostly affect z-score methods; quantile counts are rank based.
#'
#' @param data Data frame with protocol biomarkers.
#' @inheritParams validate_data
#' @param method Univariate rule: `"iqr"` (Tukey fences) or `"mad"`.
#' @param k Fence multiplier: IQRs beyond the quartiles for `"iqr"` (default 3,
#'   "far out") or robust z cut-off for `"mad"` (default 3.5).
#' @param multivariate Whether to compute robust Mahalanobis distances
#'   (requires MASS and at least `2 * p + 1` complete panels).
#' @param alpha Tail probability of the chi-square cut-off.
#' @param winsorize Whether to return a winsorized copy of `data`.
#' @param reference Optional fitted reference to measure the impact of
#'   winsorizing on scores (requires `winsorize = TRUE`).
#' @return An `al_outliers` object with `flags` (one row per measured value:
#'   fences, robust z and outlier direction), `summary` (counts per biomarker),
#'   `multivariate` (distance, cut-off and flag per complete panel), optional
#'   `winsorized` data (with attribute `al_winsorized`) and `impact` (score
#'   changes), and editable plots `univariate` and `multivariate`.
#' @family data quality
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
detect_outliers <- function(data, spec, method = c("iqr", "mad"), k = NULL, multivariate = TRUE,
                            alpha = 0.001, winsorize = FALSE, reference = NULL) {
  validate_data(data, spec); method <- match.arg(method)
  if (is.null(k)) k <- if (method == "iqr") 3 else 3.5
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) stop("k must be a positive number.", call. = FALSE)
  .check_open_unit(alpha, "alpha"); .check_flag(multivariate, "multivariate"); .check_flag(winsorize, "winsorize")
  if (!is.null(reference) && !winsorize) stop("Impact on scores requires winsorize = TRUE.", call. = FALSE)
  key <- .stratum_key(data, spec$strata)
  info <- .biomarker_info(spec)
  flags <- list(); winsorized <- data; clipped <- 0L
  for (b in spec$biomarkers) {
    x <- data[[b$name]]; tx <- .transform(x, b$transform)
    lower <- upper <- rz <- rep(NA_real_, length(x))
    for (s in unique(key)) {
      r <- which(key == s & !is.na(tx)); if (length(r) < 4L) next
      v <- tx[r]; med <- stats::median(v); mad <- stats::mad(v)
      if (method == "iqr") {
        q <- stats::quantile(v, c(0.25, 0.75), names = FALSE); iqr <- q[2] - q[1]
        lower[r] <- q[1] - k * iqr; upper[r] <- q[2] + k * iqr
      } else { lower[r] <- med - k * mad; upper[r] <- med + k * mad }
      rz[r] <- if (mad > 0) (v - med) / mad else 0
    }
    direction <- ifelse(is.na(tx) | is.na(lower), NA_character_, ifelse(tx < lower, "low", ifelse(tx > upper, "high", "none")))
    flags[[length(flags) + 1L]] <- data.frame(row = seq_along(x), biomarker = b$name, value = x,
      lower_fence = .untransform(lower, b$transform), upper_fence = .untransform(upper, b$transform),
      robust_z = rz, outlier = !is.na(direction) & direction != "none", direction = direction)
    if (winsorize) {
      out <- !is.na(direction) & direction != "none"; clipped <- clipped + sum(out)
      winsorized[[b$name]][out & direction == "low"] <- .untransform(lower, b$transform)[out & direction == "low"]
      winsorized[[b$name]][out & direction == "high"] <- .untransform(upper, b$transform)[out & direction == "high"]
    }
  }
  flags <- do.call(rbind, flags)
  summary <- .per_biomarker(spec, function(b) {
    name <- b$name
    f <- flags[flags$biomarker == name & !is.na(flags$direction), ]
    data.frame(biomarker = name, n_measured = nrow(f), n_low = sum(f$direction == "low"), n_high = sum(f$direction == "high"),
      share = if (nrow(f)) mean(f$outlier) else NA_real_)
  })
  mv <- NULL
  if (multivariate) {
    .require_suggested("MASS", "compute robust Mahalanobis distances")
    X <- vapply(spec$biomarkers, function(b) .transform(data[[b$name]], b$transform), numeric(nrow(data)))
    X <- matrix(X, nrow = nrow(data)); ok <- stats::complete.cases(X)
    if (ncol(X) >= 2L && sum(ok) >= 2L * ncol(X) + 1L) {
      fit <- .with_seed(1L, function() MASS::cov.rob(X[ok, , drop = FALSE], method = "mcd"))
      d2 <- stats::mahalanobis(X[ok, , drop = FALSE], fit$center, fit$cov)
      cutoff <- stats::qchisq(1 - alpha, df = ncol(X))
      mv <- data.frame(row = which(ok), distance = d2, cutoff = cutoff, outlier = d2 > cutoff)
    }
  }
  impact <- NULL
  if (winsorize) attr(winsorized, "al_winsorized") <- list(method = method, k = k, n_values = clipped)
  if (!is.null(reference)) {
    before <- score_load(data, reference)$scores$score; after <- score_load(winsorized, reference)$scores$score
    impact <- data.frame(row = seq_along(before), score_original = before, score_winsorized = after,
      change = after - before)
  }
  shown <- flags[!is.na(flags$value), ]
  shown$label <- factor(info$label[match(shown$biomarker, info$name)], levels = info$label)
  shown$class <- ifelse(shown$outlier, paste("Outlier", shown$direction), "Within fences")
  univariate <- ggplot2::ggplot(shown, ggplot2::aes(x = "", y = .data$value)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$class), position = ggplot2::position_jitter(width = 0.25, height = 0, seed = 3),
      size = 0.9, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(`Within fences` = .al_palette$light_grey, `Outlier high` = .al_palette$diverging[["high"]], `Outlier low` = .al_palette$categorical[[1]])) +
    ggplot2::facet_wrap(~label, scales = "free_y") + theme_allostatic() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Biomarker outliers", x = NULL, y = "Value (declared unit)", colour = NULL,
      caption = .wrap_caption(sprintf("%s rule with k = %g on the protocol transform scale%s. Values are flagged, not removed.",
        if (method == "iqr") "Tukey fence" else "Robust z (median/MAD)", k,
        if (length(spec$strata)) paste0(", within ", paste(spec$strata, collapse = ", ")) else ""), 100))
  plots <- list(univariate = univariate)
  if (!is.null(mv)) {
    mvp <- mv[order(mv$distance), ]; mvp$position <- seq_len(nrow(mvp))
    plots$multivariate <- ggplot2::ggplot(mvp, ggplot2::aes(x = .data$position, y = .data$distance, colour = .data$outlier)) +
      ggplot2::geom_hline(yintercept = mvp$cutoff[1], linetype = "22", colour = .al_palette$secondary) +
      ggplot2::geom_point(size = 1.2) + ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(values = c(`FALSE` = .al_palette$light_grey, `TRUE` = .al_palette$diverging[["high"]]), labels = c(`FALSE` = "Typical", `TRUE` = "Multivariate outlier")) +
      ggplot2::scale_x_continuous(breaks = NULL) + theme_allostatic() +
      ggplot2::labs(title = "Robust Mahalanobis distance of each complete panel", x = "Participants sorted by distance",
        y = "Squared distance (log scale)", colour = NULL,
        caption = .wrap_caption(sprintf("MCD estimate; dashed line: chi-square %g quantile with %d df.", 1 - alpha, length(info$name))))
  }
  structure(list(flags = flags, summary = summary, multivariate = mv, winsorized = if (winsorize) winsorized else NULL,
    impact = impact, method = method, k = k, plots = plots), class = "al_outliers")
}

#' @export
print.al_outliers <- function(x, ...) {
  cat(sprintf("Outlier screen (%s, k = %g)\n\n", x$method, x$k))
  print(x$summary, row.names = FALSE, digits = 3, ...)
  if (!is.null(x$multivariate)) cat(sprintf("\nMultivariate outliers: %d of %d complete panels\n", sum(x$multivariate$outlier), nrow(x$multivariate)))
  if (!is.null(x$impact)) cat(sprintf("Scores changed by winsorizing: %d\n", sum(x$impact$change != 0, na.rm = TRUE)))
  invisible(x)
}
