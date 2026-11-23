#' Build a grid of allostatic load specifications
#'
#' Generates the common operationalizations of allostatic load from one
#' protocol: contribution method, risk quantile, system balance, missing-data
#' policy and reference type (pooled, stratified or age-conditional). Invalid
#' combinations are dropped and reported, never silently altered.
#'
#' @param spec Base protocol from define_protocol(); its biomarkers and other
#'   settings (inclusive rule, bilateral rule, medications, cutoffs) are kept.
#' @param method Contribution methods to include.
#' @param risk_quantile Risk quantiles for `quantile_count` specifications.
#' @param balance Aggregation balances to include.
#' @param missing Missing-data policies to include.
#' @param min_coverage Minimum coverage used by proportion policies.
#' @param reference Reference types: `"pooled"`, `"stratified"` (by `strata`)
#'   and `"age"` (age-conditional, within `strata` when given).
#' @param strata,age Columns used by stratified and age-conditional references.
#' @return An `al_multiverse_grid` with `choices` (one row per specification),
#'   `specs` (named protocol objects) and `dropped` (invalid combinations and
#'   the reason).
#' @family multiverse analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
multiverse_grid <- function(spec, method = c("quantile_count", "z_mean"), risk_quantile = c(0.75, 0.9),
                            balance = c("biomarker", "system"), missing = "complete", min_coverage = 0.75,
                            reference = "pooled", strata = character(), age = character()) {
  validate_protocol(spec)
  method <- match.arg(method, c("quantile_count", "cutoff_count", "z_mean"), several.ok = TRUE)
  balance <- match.arg(balance, c("biomarker", "system"), several.ok = TRUE)
  missing <- match.arg(missing, c("complete", "proportion"), several.ok = TRUE)
  reference <- match.arg(reference, c("pooled", "stratified", "age"), several.ok = TRUE)
  if (!is.numeric(risk_quantile) || !length(risk_quantile) || any(!is.finite(risk_quantile)) ||
      any(risk_quantile <= 0.5 | risk_quantile >= 1)) stop("risk_quantile values must be between 0.5 and 1.", call. = FALSE)
  if ("stratified" %in% reference && !length(strata)) stop("Stratified references need strata columns.", call. = FALSE)
  if ("age" %in% reference && !length(age)) stop("Age-conditional references need an age column.", call. = FALSE)
  grid <- expand.grid(method = method, risk_quantile = sort(unique(risk_quantile)), balance = balance,
    missing = missing, reference = reference, stringsAsFactors = FALSE)
  grid$risk_quantile[grid$method != "quantile_count"] <- NA_real_
  grid <- unique(grid); rownames(grid) <- NULL
  specs <- list(); kept <- logical(nrow(grid)); reasons <- character(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    args <- utils::modifyList(unclass(spec), list(method = g$method, balance = g$balance, missing = g$missing,
      min_coverage = if (g$missing == "proportion") min_coverage else 1,
      risk_quantile = if (is.na(g$risk_quantile)) spec$risk_quantile else g$risk_quantile,
      strata = if (g$reference == "pooled") character() else strata,
      age = if (g$reference == "age") age else character()))
    built <- tryCatch(do.call(define_protocol, args), error = function(e) conditionMessage(e))
    if (is.character(built)) reasons[i] <- built else { kept[i] <- TRUE; specs[[length(specs) + 1L]] <- built }
  }
  dropped <- cbind(grid[!kept, , drop = FALSE], reason = reasons[!kept])
  grid <- grid[kept, , drop = FALSE]
  grid$spec_id <- sprintf("S%02d", seq_len(nrow(grid)))
  grid$label <- paste(grid$method, ifelse(is.na(grid$risk_quantile), "", paste0("q", grid$risk_quantile)),
    grid$balance, grid$missing, grid$reference)
  grid$label <- gsub(" +", " ", grid$label)
  names(specs) <- grid$spec_id
  rownames(grid) <- NULL; rownames(dropped) <- NULL
  structure(list(choices = grid[c("spec_id", "method", "risk_quantile", "balance", "missing", "reference", "label")],
    specs = specs, dropped = dropped), class = "al_multiverse_grid")
}

.kappa <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  if (length(a) < 2L) return(NA_real_)
  po <- mean(a == b); pe <- mean(a) * mean(b) + mean(!a) * mean(!b)
  if (pe >= 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

#' Multiverse analysis of allostatic load specifications
#'
#' Scores the same participants under every specification of a grid, then
#' describes how robust each participant's classification is, how strongly
#' specifications agree and, optionally, how an outcome association varies
#' across specifications (a specification curve). It describes robustness to
#' analytic choices; selecting a specification by its outcome association
#' would bias inference.
#'
#' @param data Study data frame.
#' @param reference_data Calibration data for every specification (may be the
#'   study data itself, which then calibrates each specification internally).
#' @param grid An `al_multiverse_grid` from multiverse_grid().
#' @param id,group Optional identifier and group columns.
#' @param outcome Optional outcome column for the specification curve.
#' @param covariates Optional adjustment columns for the outcome models.
#' @param family `"gaussian"` (linear model) or `"binomial"` (logistic model).
#' @param top Fraction of highest scores defining "high load" in each
#'   specification. Percentiles use mid-ranks, so tied scores share one
#'   classification: a tie group is high when its mid-rank percentile
#'   exceeds `1 - top`. Coarse count scores therefore flag about, not
#'   exactly, the requested fraction.
#' @return An `al_multiverse` object with `choices`, `results` (one al_result
#'   per specification), `failures`, `scores` (long format with within-
#'   specification percentiles and high-load flags), `stability` (per
#'   participant: share of specifications classifying them as high load and
#'   class `robustly high`, `robustly not high` or `specification-dependent`,
#'   using 90% and 10% limits), `agreement` (pairwise Spearman correlation,
#'   kappa and share classified alike), `effects` (outcome association per SD
#'   of each score, with 95% Wald intervals) and editable plots.
#' @references Steegen S et al. (2016). Increasing transparency through a
#'   multiverse analysis. Perspectives on Psychological Science 11, 702-712.
#'   McLoughlin S, Kenny RA, McCrory C (2020). Does the choice of allostatic
#'   load scoring algorithm matter for predicting age-related health outcomes?
#'   Psychoneuroendocrinology 120, 104789.
#' @family multiverse analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
multiverse_analysis <- function(data, reference_data, grid, id = NULL, group = NULL, outcome = NULL,
                                covariates = character(), family = c("gaussian", "binomial"), top = 0.25) {
  if (!inherits(grid, "al_multiverse_grid")) stop("grid must come from multiverse_grid().", call. = FALSE)
  if (length(grid$specs) < 2L) stop("The grid must contain at least two specifications.", call. = FALSE)
  family <- match.arg(family); .check_open_unit(top, "top")
  if (!is.null(outcome)) {
    .scalar_string(outcome, "outcome")
    if (!all(c(outcome, covariates) %in% names(data))) stop("Missing outcome or covariate columns.", call. = FALSE)
  }
  results <- list(); failures <- data.frame(spec_id = character(), error = character())
  for (k in names(grid$specs)) {
    r <- tryCatch(score_load(data, fit_reference(reference_data, grid$specs[[k]], label = k), id = id, group = group),
      error = function(e) conditionMessage(e))
    if (is.character(r)) failures <- rbind(failures, data.frame(spec_id = k, error = r)) else results[[k]] <- r
  }
  if (length(results) < 2L) stop("Fewer than two specifications could be scored: ", paste(failures$error, collapse = "; "), call. = FALSE)
  choices <- grid$choices[grid$choices$spec_id %in% names(results), ]
  scores <- do.call(rbind, lapply(names(results), function(k) {
    s <- results[[k]]$scores; x <- s$score; valid <- !is.na(x)
    percentile <- rep(NA_real_, length(x)); percentile[valid] <- (rank(x[valid]) - 0.5) / sum(valid)
    data.frame(spec_id = k, row = s$row, id = s$id, score = x, percentile = percentile,
      high = ifelse(valid, percentile > 1 - top, NA))
  }))
  high <- matrix(scores$high, ncol = length(results)); pct <- matrix(scores$percentile, ncol = length(results))
  n_valid <- rowSums(!is.na(high))
  share_high <- ifelse(n_valid > 0, rowSums(high, na.rm = TRUE) / pmax(n_valid, 1), NA_real_)
  stability <- data.frame(row = results[[1]]$scores$row, id = results[[1]]$scores$id, n_valid_specs = n_valid,
    share_high = share_high,
    median_percentile = apply(pct, 1, function(p) if (all(is.na(p))) NA_real_ else stats::median(p, na.rm = TRUE)),
    min_percentile = apply(pct, 1, function(p) if (all(is.na(p))) NA_real_ else min(p, na.rm = TRUE)),
    max_percentile = apply(pct, 1, function(p) if (all(is.na(p))) NA_real_ else max(p, na.rm = TRUE)),
    class = ifelse(is.na(share_high), NA_character_, ifelse(share_high >= 0.9, "robustly high",
      ifelse(share_high <= 0.1, "robustly not high", "specification-dependent"))))
  if (!is.null(group)) stability$group <- results[[1]]$scores$group
  score_mat <- matrix(scores$score, ncol = length(results), dimnames = list(NULL, names(results)))
  pairs <- utils::combn(length(results), 2)
  agreement <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    a <- pairs[1, j]; b <- pairs[2, j]; ok <- !is.na(score_mat[, a]) & !is.na(score_mat[, b])
    data.frame(spec_a = names(results)[a], spec_b = names(results)[b], n_matched = sum(ok),
      spearman = .safe_spearman(score_mat[, a], score_mat[, b]), kappa = .kappa(high[, a], high[, b]),
      share_same_class = if (any(ok)) mean(high[ok, a] == high[ok, b]) else NA_real_)
  }))
  effects <- NULL
  if (!is.null(outcome)) effects <- do.call(rbind, lapply(names(results), function(k) {
    d <- data.frame(y = data[[outcome]], score = results[[k]]$scores$score, data[covariates])
    d <- d[stats::complete.cases(d), , drop = FALSE]
    if (nrow(d) < length(covariates) + 3L || stats::sd(d$score) == 0)
      return(data.frame(spec_id = k, n = nrow(d), estimate = NA_real_, std_error = NA_real_, lower = NA_real_, upper = NA_real_))
    d$score <- (d$score - mean(d$score)) / stats::sd(d$score)
    formula <- stats::reformulate(c("score", if (length(covariates)) paste0("`", covariates, "`")), response = "y")
    fit <- stats::glm(formula, data = d, family = family)
    est <- stats::coef(summary(fit))["score", ]
    data.frame(spec_id = k, n = nrow(d), estimate = est[[1]], std_error = est[[2]],
      lower = est[[1]] - stats::qnorm(0.975) * est[[2]], upper = est[[1]] + stats::qnorm(0.975) * est[[2]])
  }))
  if (!is.null(effects)) effects <- merge(choices, effects, by = "spec_id", sort = FALSE)
  classes <- table(factor(stability$class, c("robustly high", "specification-dependent", "robustly not high")))
  summary <- data.frame(n_specs = length(results), n_failed = nrow(failures), n_participants = nrow(stability),
    robustly_high = classes[[1]], specification_dependent = classes[[2]], robustly_not_high = classes[[3]],
    median_spearman = stats::median(agreement$spearman, na.rm = TRUE), median_kappa = stats::median(agreement$kappa, na.rm = TRUE),
    top = top)
  if (!is.null(effects)) {
    summary$median_estimate <- stats::median(effects$estimate, na.rm = TRUE)
    summary$share_ci_excludes_zero <- mean(effects$lower > 0 | effects$upper < 0, na.rm = TRUE)
  }
  out <- structure(list(choices = choices, results = results, failures = failures, scores = scores,
    stability = stability, agreement = agreement, effects = effects, summary = summary, family = family,
    outcome = outcome, plots = list()), class = "al_multiverse")
  out$plots <- .multiverse_plots(out)
  out
}

.multiverse_plots <- function(x) {
  plots <- list()
  class_colors <- c(`robustly high` = .al_palette$diverging[["high"]], `specification-dependent` = .al_palette$categorical[[4]], `robustly not high` = .al_palette$light_grey)
  st <- x$stability[!is.na(x$stability$class), ]
  note <- ""
  if (nrow(st) > 500L) {
    st <- st[.with_seed(1L, function() sort(sample.int(nrow(st), 500L))), ]
    note <- sprintf(" Random sample of 500 of %d participants.", sum(!is.na(x$stability$class)))
  }
  if (nrow(st)) {
    st <- st[order(st$median_percentile), ]; st$position <- seq_len(nrow(st))
    counts <- table(factor(x$stability$class, names(class_colors)))
    labels <- paste0(names(class_colors), " (n=", as.integer(counts), ")")
    plots$stability <- .finish_plot(
      ggplot2::ggplot(st, ggplot2::aes(x = .data$position)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = .data$min_percentile, ymax = .data$max_percentile, colour = .data$class),
          linewidth = 0.5) +
        ggplot2::geom_point(ggplot2::aes(y = .data$median_percentile), size = 0.7, colour = .al_palette$primary) +
        ggplot2::geom_hline(yintercept = 1 - x$summary$top, colour = .al_palette$secondary, linetype = "22", linewidth = 0.4) +
        ggplot2::scale_colour_manual(values = class_colors, labels = labels, breaks = names(class_colors)) +
        ggplot2::scale_y_continuous(labels = .percent, limits = c(0, 1)) + ggplot2::scale_x_continuous(breaks = NULL),
      title = "Who is high-load under every specification?",
      subtitle = sprintf("Percentile range across %d specifications; dashed: top-%s cut",
        x$summary$n_specs, .percent(x$summary$top)),
      x = "Participants sorted by median percentile", y = "Percentile of score", colour = NULL, width = 100,
      caption = paste0("Robust: classified alike in at least 90% of valid specifications.", note))
  }
  ids <- x$choices$spec_id
  m <- rbind(x$agreement[c("spec_a", "spec_b", "spearman")],
    stats::setNames(x$agreement[c("spec_b", "spec_a", "spearman")], c("spec_a", "spec_b", "spearman")),
    data.frame(spec_a = ids, spec_b = ids, spearman = 1))
  m$spec_a <- factor(m$spec_a, ids); m$spec_b <- factor(m$spec_b, rev(ids))
  m$label_text <- sprintf("%.2f", m$spearman)
  heat <- .heatmap_layer(m, "spec_a", "spec_b", "spearman", label = "label_text", threshold = 0.75,
    label_size = 2.6, colour = "white", linewidth = 0.5)
  agreement <- .finish_plot(
    ggplot2::ggplot(m, ggplot2::aes(x = .data$spec_a, y = .data$spec_b, fill = .data$spearman)) +
      heat$tile +
      ggplot2::scale_fill_gradient(low = .al_palette$light_blue, high = .al_palette$categorical[[1]], limits = c(min(0, m$spearman, na.rm = TRUE), 1),
        na.value = .al_palette$unavailable),
    title = "Agreement between specifications", x = NULL, y = NULL, fill = "Spearman", width = 110,
    caption = paste(paste0(x$choices$spec_id, ": ", x$choices$label), collapse = "; ")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (length(ids) <= 16L) agreement <- agreement + heat$text
  plots$agreement <- agreement
  if (!is.null(x$effects) && any(is.finite(x$effects$estimate))) {
    e <- x$effects[is.finite(x$effects$estimate), ]
    e <- e[order(e$estimate), ]; e$position <- seq_len(nrow(e))
    e$excludes <- ifelse(e$lower > 0 | e$upper < 0, "95% CI excludes 0", "95% CI includes 0")
    curve <- .finish_plot(
      ggplot2::ggplot(e, ggplot2::aes(x = .data$position, y = .data$estimate)) +
        ggplot2::geom_hline(yintercept = 0, colour = .al_palette$secondary, linewidth = 0.4) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), colour = .al_palette$light_grey, linewidth = 0.7) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$excludes), size = 2) +
        ggplot2::scale_colour_manual(values = c(`95% CI excludes 0` = .al_palette$categorical[[1]], `95% CI includes 0` = .al_palette$muted)) +
        ggplot2::scale_x_continuous(breaks = NULL, expand = ggplot2::expansion(add = 0.6)),
      title = "Specification curve", x = NULL, colour = NULL,
      subtitle = sprintf("Association of %s with each score (per SD; %s model)", x$outcome,
        if (x$family == "binomial") "logistic, log-odds" else "linear"),
      y = "Estimate per SD of score")
    dims <- c("method", "risk_quantile", "balance", "missing", "reference")
    ind <- do.call(rbind, lapply(dims, function(d) data.frame(position = e$position, dimension = d,
      level = ifelse(is.na(e[[d]]), "not applicable", as.character(e[[d]])))))
    ind <- ind[ind$level != "not applicable", ]
    ind$row <- paste(ind$dimension, ind$level, sep = ": ")
    ind$row <- factor(ind$row, levels = rev(unique(ind$row[order(match(ind$dimension, dims), ind$level)])))
    choice_plot <- .finish_plot(
      ggplot2::ggplot(ind, ggplot2::aes(x = .data$position, y = .data$row)) +
        ggplot2::geom_point(size = 1.8, colour = .al_palette$secondary) +
        ggplot2::scale_x_continuous(breaks = NULL, limits = range(e$position), expand = ggplot2::expansion(add = 0.6)),
      x = "Specifications sorted by estimate", y = NULL)
    plots$curve <- curve; plots$choices <- choice_plot
    if (requireNamespace("patchwork", quietly = TRUE))
      plots$specification_curve <- patchwork::wrap_plots(curve, choice_plot, ncol = 1, heights = c(1.4, 1))
  }
  plots
}

#' @export
print.al_multiverse <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Allostatic load multiverse: %d specifications (%d failed), %d participants\n",
    s$n_specs, s$n_failed, s$n_participants))
  cat(sprintf("Top %s classification: robustly high %d | specification-dependent %d | robustly not high %d\n",
    .percent(s$top), s$robustly_high, s$specification_dependent, s$robustly_not_high))
  cat(sprintf("Median Spearman %.2f | median kappa %.2f\n", s$median_spearman, s$median_kappa))
  if (!is.null(x$effects)) cat(sprintf("Outcome %s: median estimate per SD %.3f; 95%% CI excludes 0 in %s of specifications\n",
    x$outcome, s$median_estimate, .percent(s$share_ci_excludes_zero)))
  invisible(x)
}
