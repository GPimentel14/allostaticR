.check_variation <- function(variation, spec) {
  names_all <- .biomarker_field(spec, "name")
  needed <- c("biomarker", "cv_analytical", "cv_within", "source")
  if (!is.data.frame(variation) || !all(needed %in% names(variation)))
    stop("variation must be a data frame with columns biomarker, cv_analytical, cv_within and source.", call. = FALSE)
  if (!is.character(variation$biomarker) || anyDuplicated(variation$biomarker) || !setequal(variation$biomarker, names_all))
    stop("variation must list every protocol biomarker exactly once.", call. = FALSE)
  for (column in c("cv_analytical", "cv_within")) {
    x <- variation[[column]]
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0) || any(x >= 1))
      stop(column, " must contain coefficients of variation as proportions in [0, 1).", call. = FALSE)
  }
  if (!is.character(variation$source) || anyNA(variation$source) || any(!nzchar(trimws(variation$source))))
    stop("Every coefficient of variation needs a non-empty source.", call. = FALSE)
  out <- variation[match(names_all, variation$biomarker), needed]
  out$sigma <- sqrt(log(1 + out$cv_analytical^2 + out$cv_within^2))
  rownames(out) <- NULL
  out
}

# Probability that a log-normal quantity lies beyond a threshold.
.beyond <- function(m, s, threshold, upper) {
  out <- rep(0, length(m)); ok <- is.finite(threshold)
  if (!any(ok)) return(out)
  lt <- log(threshold[ok])
  out[ok] <- ifelse(s[ok] > 0, if (upper) stats::pnorm(lt, m[ok], s[ok], lower.tail = FALSE) else stats::pnorm(lt, m[ok], s[ok]),
    if (upper) as.numeric(m[ok] >= lt) else as.numeric(m[ok] <= lt))
  out
}

# Exact distribution of a sum of independent Bernoulli variables (Poisson-binomial).
.poisson_binomial <- function(p) {
  d <- 1
  for (q in p) d <- c(d * (1 - q), 0) + c(0, d * q)
  d
}

# Exact distribution of a score given component probabilities; returns list(value, prob).
.score_distribution <- function(p, systems, scale) {
  if (scale %in% c("count", "proportion")) {
    prob <- .poisson_binomial(p); k <- seq_along(prob) - 1L
    return(list(value = if (scale == "count") k else k / length(p), prob = prob))
  }
  value <- 0; prob <- 1; levels <- unique(systems)
  for (s in levels) {
    q <- p[systems == s]; pb <- .poisson_binomial(q); v <- (seq_along(pb) - 1L) / length(q) / length(levels)
    value <- as.vector(outer(value, v, `+`)); prob <- as.vector(outer(prob, pb))
    key <- round(value, 10); prob <- tapply(prob, key, sum); value <- as.numeric(names(prob)); prob <- as.vector(prob)
  }
  list(value = value, prob = prob)
}

.distribution_quantile <- function(d, q) d$value[order(d$value)][which(cumsum(d$prob[order(d$value)]) >= q - 1e-12)[1]]

#' Uncertainty from analytical and within-person biological variation
#'
#' A single measurement close to a threshold is a near coin-flip: with
#' analytical and within-person biological variation, the person's
#' homeostatic set point may lie on either side. Given coefficients of
#' variation for every biomarker, this function returns the probability that
#' each biomarker is truly at risk and the exact distribution of each
#' participant's load score.
#'
#' Measurements are modelled as log-normal around the set point with
#' log-scale SD `sqrt(log(1 + CVa^2 + CVi^2))`. With `prior = "flat"`, the
#' set point is centred on the measurement (the metrological convention for
#' decision limits). With `prior = "empirical"`, the observed between-person
#' distribution of the study sample shrinks the set point towards the
#' population mean, which also captures regression to the mean; it requires
#' the between-person variance to exceed the within-person variance.
#' Biomarkers are treated as independent given their set points. A TRUE
#' medication flag makes a biomarker certainly at risk.
#'
#' @param result An `al_result` from a count protocol.
#' @param variation Data frame with one row per protocol biomarker and columns
#'   `biomarker`, `cv_analytical`, `cv_within` (coefficients of variation as
#'   proportions, for example 0.05 for 5%) and `source` (the citation of each
#'   value, for example the EFLM Biological Variation Database). No default
#'   values are supplied by the package.
#' @param conf Coverage of the per-participant score interval.
#' @param prior `"flat"` or `"empirical"`; see Details.
#' @param borderline Probability range regarded as borderline classification.
#' @return An `al_uncertainty` object with `components` (probability that each
#'   measured biomarker is truly at risk), `scores` (observed score, expected
#'   score, interval and probability that the true score equals the observed
#'   one), `biomarkers` (share of borderline classifications), `summary`,
#'   `variation` (with the log-scale SD) and editable plots: `intervals`,
#'   `borderline` and, for counts, `reclassification` (expected distribution of
#'   true scores for each observed score). These describe measurement and
#'   within-person variation conditional on fixed reference thresholds; for
#'   reference-calibration uncertainty see [bootstrap_load()].
#' @references Fraser CG (2001). Biological Variation: From Principles to
#'   Practice. AACC Press. Aarsand AK et al. (2020). The EFLM Biological
#'   Variation Database. https://biologicalvariation.eu/
#' @examples
#' source(system.file("examples", "panel.R", package = "allostaticR"))
#' result <- allostatic_load(data, spec)
#' # Illustrative values only: supply sourced CVs for real analyses.
#' cv <- data.frame(biomarker = c("sbp", "hdl", "crp"), cv_analytical = c(0.02, 0.02, 0.05),
#'   cv_within = c(0.08, 0.07, 0.40), source = "Illustrative values, not from a database")
#' u <- measurement_uncertainty(result, cv)
#' u$summary
#' @family uncertainty
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
measurement_uncertainty <- function(result, variation, conf = 0.90, prior = c("flat", "empirical"),
                                    borderline = c(0.1, 0.9)) {
  .assert_result(result); prior <- match.arg(prior); .check_open_unit(conf, "conf")
  spec <- result$reference$spec
  if (!.is_count_method(spec)) stop("Measurement uncertainty requires a count protocol.", call. = FALSE)
  if (!is.numeric(borderline) || length(borderline) != 2L || any(!is.finite(borderline)) ||
      borderline[1] <= 0 || borderline[2] >= 1 || borderline[1] >= borderline[2])
    stop("borderline must be two increasing probabilities in (0, 1).", call. = FALSE)
  variation <- .check_variation(variation, spec)
  cmp <- result$components
  if (any(cmp$value <= 0, na.rm = TRUE)) stop("The log-normal variation model requires positive values.", call. = FALSE)
  sigma <- variation$sigma[match(cmp$biomarker, variation$biomarker)]
  values <- if (spec$medication_rule == "add") cmp$value_adjusted else cmp$value
  if (any(values <= 0, na.rm = TRUE)) stop("The log-normal variation model requires positive values.", call. = FALSE)
  m <- log(values); s <- sigma
  if (prior == "empirical") for (b in variation$biomarker) {
    i <- cmp$biomarker == b & is.finite(m); lx <- m[i]
    tau2 <- stats::var(lx) - variation$sigma[variation$biomarker == b]^2
    if (!is.finite(tau2) || tau2 <= 0)
      stop("Within-person variation exceeds the observed between-person variance for ", b, "; use prior = 'flat'.", call. = FALSE)
    k <- tau2 / (tau2 + variation$sigma[variation$biomarker == b]^2)
    m[i] <- mean(lx) + k * (lx - mean(lx)); s[i] <- sqrt(k) * variation$sigma[variation$biomarker == b]
  }
  p <- .beyond(m, s, cmp$lower, upper = FALSE) * (cmp$direction %in% c("low", "both")) +
    .beyond(m, s, cmp$upper, upper = TRUE) * (cmp$direction %in% c("high", "both"))
  if (spec$medication_rule == "flag") p[cmp$medication_flag] <- 1
  p[is.na(cmp$risk)] <- NA_real_
  components <- data.frame(row = cmp$row, id = cmp$id, biomarker = cmp$biomarker, system = cmp$system,
    value = cmp$value, lower = cmp$lower, upper = cmp$upper, risk = cmp$risk, p_risk = p,
    borderline = !is.na(p) & p >= borderline[1] & p <= borderline[2])
  n <- nrow(result$scores); pm <- matrix(p, nrow = n)
  systems <- .biomarker_field(spec, "system")
  scale <- result$audit$scale
  rows <- which(!is.na(result$scores$score))
  per_person <- lapply(rows, function(i) {
    observed <- !is.na(pm[i, ])
    d <- .score_distribution(pm[i, observed], systems[observed], scale)
    c(expected = sum(d$value * d$prob), lower = .distribution_quantile(d, (1 - conf) / 2),
      upper = .distribution_quantile(d, 1 - (1 - conf) / 2),
      p_observed = sum(d$prob[abs(d$value - result$scores$score[i]) < 1e-9]))
  })
  stats_mat <- if (length(per_person)) do.call(rbind, per_person) else matrix(numeric(), 0, 4,
    dimnames = list(NULL, c("expected", "lower", "upper", "p_observed")))
  scores <- data.frame(row = result$scores$row, id = result$scores$id, score = result$scores$score,
    expected = NA_real_, lower = NA_real_, upper = NA_real_, p_observed = NA_real_,
    n_borderline = rowSums(matrix(components$borderline, nrow = n)))
  if ("group" %in% names(result$scores)) scores$group <- result$scores$group
  scores[rows, c("expected", "lower", "upper", "p_observed")] <- stats_mat
  biomarkers <- .per_biomarker(spec, function(b) {
    name <- b$name
    d <- components[components$biomarker == name & !is.na(components$p_risk), ]
    data.frame(biomarker = name, system = systems[match(name, variation$biomarker)], n_measured = nrow(d),
      share_borderline = if (nrow(d)) mean(d$borderline) else NA_real_,
      mean_p_flagged = .mean_or_na(d$p_risk[d$risk == 1L]), mean_p_not_flagged = .mean_or_na(d$p_risk[d$risk == 0L]))
  })
  valid <- scores[rows, ]
  summary <- data.frame(n_scored = length(rows), mean_observed = .mean_or_na(valid$score),
    mean_expected = .mean_or_na(valid$expected), mean_p_observed = .mean_or_na(valid$p_observed),
    share_with_borderline = if (length(rows)) mean(valid$n_borderline > 0) else NA_real_,
    share_interval_excludes_observed = if (length(rows)) mean(valid$score < valid$lower | valid$score > valid$upper) else NA_real_,
    conf = conf, prior = prior, scale = scale)
  out <- structure(list(components = components, scores = scores, biomarkers = biomarkers, summary = summary,
    variation = variation, plots = list()), class = "al_uncertainty")
  out$plots <- .uncertainty_plots(out, spec, pm, rows)
  out
}

.uncertainty_plots <- function(x, spec, pm, rows) {
  colors <- .system_colors(.protocol_systems(spec)); s <- x$summary
  plots <- list()
  valid <- x$scores[rows, ]
  if (nrow(valid)) {
    shown <- if (nrow(valid) > 300L) valid[.with_seed(1L, function() sort(sample.int(nrow(valid), 300L))), ] else valid
    shown <- shown[order(shown$expected, shown$score), ]; shown$position <- seq_len(nrow(shown))
    plots$intervals <- ggplot2::ggplot(shown, ggplot2::aes(x = .data$position)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), colour = .al_palette$light_grey, linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(y = .data$expected, shape = "Expected"), colour = .al_palette$secondary, size = 1.2) +
      ggplot2::geom_point(ggplot2::aes(y = .data$score, shape = "Observed"), colour = .al_palette$categorical[[1]], size = 1.4) +
      ggplot2::scale_shape_manual(values = c(Expected = 16, Observed = 4)) +
      ggplot2::scale_x_continuous(breaks = NULL) + theme_allostatic() +
      ggplot2::labs(title = "Scores under measurement and biological variation",
        subtitle = sprintf("Gray: %g%% interval of the set-point score; %s", 100 * s$conf,
          if (nrow(valid) > 300L) sprintf("random sample of 300 of %d participants", nrow(valid)) else "all scored participants"),
        x = "Participants sorted by expected score", y = .scale_label(s$scale), shape = NULL,
        caption = .wrap_caption(sprintf("Log-normal model with supplied CVs (%s prior). Mean observed %.2f, mean expected %.2f; %s of participants have a borderline biomarker.",
          s$prior, s$mean_observed, s$mean_expected, .percent(s$share_with_borderline)), 100))
  }
  b <- x$biomarkers; b$label <- factor(b$biomarker, levels = rev(b$biomarker))
  plots$borderline <- ggplot2::ggplot(b, ggplot2::aes(x = .data$share_borderline, y = .data$label, fill = .data$system)) +
    ggplot2::geom_col(width = 0.65) + ggplot2::scale_fill_manual(values = colors, breaks = names(colors)) +
    ggplot2::scale_x_continuous(labels = .percent, limits = c(0, max(0.1, b$share_borderline, na.rm = TRUE)),
      expand = ggplot2::expansion(mult = c(0, 0.05))) + theme_allostatic() +
    ggplot2::labs(title = "Borderline risk classifications", x = "Measured values with uncertain classification",
      y = NULL, fill = "System", caption = "Borderline: probability of being truly at risk between the chosen limits (default 10% and 90%).")
  if (s$scale == "count" && nrow(valid)) {
    j <- ncol(pm)
    mass <- do.call(rbind, lapply(rows, function(i) data.frame(observed = x$scores$score[i], true = 0:j,
      prob = .poisson_binomial(pm[i, ]))))
    tab <- stats::aggregate(prob ~ observed + true, data = mass, FUN = sum)
    tab$share <- tab$prob / stats::ave(tab$prob, tab$observed, FUN = sum)
    tab <- tab[tab$share >= 0.005, ]
    plots$reclassification <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$observed, y = .data$true, fill = .data$share)) +
      ggplot2::geom_tile(width = 0.94, height = 0.94) +
      ggplot2::geom_text(ggplot2::aes(label = .percent(.data$share)), size = 3,
        colour = ifelse(tab$share > 0.55, "white", .al_palette$primary)) +
      ggplot2::scale_fill_gradient(low = .al_palette$light_blue, high = .al_palette$categorical[[1]], limits = c(0, 1), labels = .percent) +
      ggplot2::scale_x_continuous(breaks = 0:j) + ggplot2::scale_y_continuous(breaks = 0:j) +
      ggplot2::coord_equal() + theme_allostatic() +
      ggplot2::labs(title = "How certain is each observed load?", x = "Observed allostatic load (count)",
        y = "Set-point load (count)", fill = "Share",
        caption = .wrap_caption("Each column: expected distribution of set-point loads among participants with that observed load.", 90))
  }
  plots
}

#' @export
print.al_uncertainty <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Measurement uncertainty (%s prior, %s scale): %d scored participants\n", s$prior, s$scale, s$n_scored))
  cat(sprintf("Mean observed %.3f | mean expected %.3f | mean P(true = observed) %.2f\n",
    s$mean_observed, s$mean_expected, s$mean_p_observed))
  cat(sprintf("Participants with a borderline biomarker: %s\n\n", .percent(s$share_with_borderline)))
  print(x$biomarkers, row.names = FALSE, digits = 3, ...)
  invisible(x)
}

.rcv_changes <- function(trajectory, out, variation, conf) {
  spec <- trajectory$reference$spec
  if (!.is_count_method(spec)) stop("RCV classification of risk changes requires a count protocol.", call. = FALSE)
  .check_open_unit(conf, "conf")
  variation <- .check_variation(variation, spec)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  cmp <- trajectory$components
  weights <- stats::setNames(.score_weights(spec, trajectory$scale), .biomarker_field(spec, "name"))
  details <- list(); out$n_flips <- out$n_confirmed_flips <- NA_integer_; out$confirmed_change <- NA_real_
  for (k in which(out$status == "ok")) {
    d <- cmp[cmp$id == out$id[k], ]
    first <- d[d$time == out$first_time[k], ]; last <- d[d$time == out$last_time[k], ]
    last <- last[match(first$biomarker, last$biomarker), ]
    sigma <- variation$sigma[match(first$biomarker, variation$biomarker)]
    log_change <- log(last$value / first$value)
    limit <- z * sqrt(2) * sigma
    flip <- !is.na(first$risk) & !is.na(last$risk) & first$risk != last$risk
    reason <- ifelse(!flip, "none", ifelse(first$medication_flag != last$medication_flag, "medication", "measurement"))
    beyond <- is.finite(log_change) & abs(log_change) > limit
    confirmed <- flip & reason == "measurement" & beyond
    complete <- all(!is.na(first$risk)) && all(!is.na(last$risk))
    out$n_flips[k] <- sum(flip); out$n_confirmed_flips[k] <- sum(confirmed)
    out$confirmed_change[k] <- if (complete) sum((last$risk - first$risk)[confirmed] * weights[first$biomarker[confirmed]]) else NA_real_
    details[[length(details) + 1L]] <- data.frame(id = out$id[k], biomarker = first$biomarker,
      first_value = first$value, last_value = last$value, first_risk = first$risk, last_risk = last$risk,
      log_change = log_change, rcv_up = exp(limit) - 1, rcv_down = exp(-limit) - 1, beyond_rcv = beyond,
      flip = flip, reason = reason, confirmed = confirmed)
  }
  attr(out, "biomarker_changes") <- if (length(details)) do.call(rbind, details) else NULL
  out
}
