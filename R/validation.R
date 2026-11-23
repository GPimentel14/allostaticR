.auc <- function(y, p) {
  ok <- !is.na(y) & is.finite(p); y <- as.logical(y[ok]); p <- p[ok]
  n1 <- sum(y); n0 <- sum(!y)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p)[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.fold_ids <- function(y, folds, type) {
  n <- length(y); id <- integer(n)
  groups <- if (type == "binary") split(seq_len(n), as.logical(y)) else list(seq_len(n))
  for (g in groups) id[g] <- sample(rep_len(seq_len(folds), length(g)))
  id
}

#' Cross-validated comparison of scoring algorithms against an outcome
#'
#' Compares allostatic load specifications by how well each score predicts an
#' outcome in held-out data, using repeated k-fold cross-validation. Unless an
#' external reference sample is supplied, reference thresholds are refitted in
#' every training fold, so test participants never inform their own
#' thresholds. This is the appropriate way to ask which algorithm predicts
#' best; selecting an algorithm on in-sample fit would be optimistic.
#'
#' @param data Data frame with biomarkers, outcome and covariates.
#' @param specs Named list of protocols (at least one), or an
#'   `al_multiverse_grid`.
#' @param outcome Outcome column: 0/1 or logical for `"binary"`, numeric for
#'   `"continuous"`, event indicator for `"survival"`.
#' @param type Outcome type.
#' @param time Follow-up time column for survival outcomes.
#' @param covariates Adjustment columns. When supplied, a covariate-only model
#'   is also validated and incremental performance is reported.
#' @param reference_data Optional external reference sample; when NULL the
#'   reference is fitted in each training fold.
#' @param folds,repeats Number of folds and repetitions.
#' @inheritParams simulate_biomarkers
#' @return An `al_validation` object with `folds` (metrics for every
#'   specification, repetition and fold), `summary` (mean and SD across folds:
#'   AUC and Brier score for binary outcomes, Harrell's C for survival, R^2 and
#'   RMSE for continuous outcomes, plus the calibration slope and the
#'   difference from the covariate-only model), `failures` (data frame of any
#'   failed fold/specification evaluations and their error messages) and an
#'   editable plot. Scores are standardized within each training fold.
#' @family multiverse analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
validate_scores <- function(data, specs, outcome, type = c("binary", "continuous", "survival"), time = NULL,
                            covariates = character(), reference_data = NULL, folds = 5L, repeats = 1L, seed = 1L) {
  type <- match.arg(type)
  if (inherits(specs, "al_multiverse_grid")) specs <- specs$specs
  .named_list(specs); invisible(lapply(specs, validate_protocol))
  .check_integer(folds, "folds", 2L, 20L); .check_integer(repeats, "repeats", 1L, 50L)
  .scalar_string(outcome, "outcome")
  needed <- c(outcome, time, covariates)
  if (!all(needed %in% names(data))) stop("Missing outcome, time or covariate columns.", call. = FALSE)
  if (type == "survival") {
    .require_suggested("survival", "validate survival outcomes")
    if (is.null(time)) stop("Survival outcomes need a time column.", call. = FALSE)
  }
  y <- data[[outcome]]
  if (type == "binary" && !all(stats::na.omit(y) %in% c(0, 1, TRUE, FALSE))) stop("Binary outcomes must be 0/1 or logical.", call. = FALSE)
  keep <- stats::complete.cases(data[needed])
  data <- data[keep, , drop = FALSE]; y <- data[[outcome]]
  if (nrow(data) < 4L * folds) stop("Too few complete rows for cross-validation.", call. = FALSE)
  fit_predict <- function(train, test, s_train, s_test) {
    standardize <- !is.null(s_train)
    if (standardize) { mu <- mean(s_train); sdv <- stats::sd(s_train); train$.score <- (s_train - mu) / sdv; test$.score <- (s_test - mu) / sdv }
    terms <- c(if (standardize) ".score", if (length(covariates)) paste0("`", covariates, "`"))
    if (!length(terms)) terms <- "1"
    if (type == "survival") {
      f <- stats::as.formula(paste0("survival::Surv(`", time, "`, `", outcome, "`) ~ ", paste(terms, collapse = " + ")))
      model <- survival::coxph(f, data = train)
      lp <- if (length(stats::coef(model))) stats::predict(model, newdata = test, type = "lp") else rep(0, nrow(test))
      c_index <- survival::concordance(stats::as.formula(paste0("survival::Surv(`", time, "`, `", outcome, "`) ~ lp")),
        data = cbind(test, lp = lp), reverse = TRUE)$concordance
      slope <- tryCatch(unname(stats::coef(survival::coxph(stats::as.formula(paste0("survival::Surv(`", time, "`, `", outcome, "`) ~ lp")),
        data = cbind(test, lp = lp)))), error = function(e) NA_real_)
      return(c(c_index = c_index, calibration_slope = slope))
    }
    f <- stats::as.formula(paste0("`", outcome, "` ~ ", paste(terms, collapse = " + ")))
    model <- stats::glm(f, data = train, family = if (type == "binary") stats::binomial() else stats::gaussian())
    lp <- stats::predict(model, newdata = test, type = "link")
    obs <- as.numeric(test[[outcome]])
    slope <- tryCatch(unname(stats::coef(stats::glm(obs ~ lp, family = if (type == "binary") stats::binomial() else stats::gaussian()))[2]),
      error = function(e) NA_real_)
    if (type == "binary") { p <- stats::plogis(lp); return(c(auc = .auc(obs, p), brier = mean((p - obs)^2), calibration_slope = slope)) }
    c(r_squared = 1 - sum((obs - lp)^2) / sum((obs - mean(obs))^2), rmse = sqrt(mean((obs - lp)^2)), calibration_slope = slope)
  }
  failures <- list()
  rows <- .with_seed(seed, function() {
    out <- list()
    for (r in seq_len(repeats)) {
      fold <- .fold_ids(if (type == "survival") data[[outcome]] else y, folds, type)
      for (k in seq_len(folds)) {
        train <- data[fold != k, , drop = FALSE]; test <- data[fold == k, , drop = FALSE]
        if (length(covariates)) {
          m <- tryCatch(fit_predict(train, test, NULL, NULL), error = function(e) NULL)
          if (!is.null(m)) out[[length(out) + 1L]] <- data.frame(spec = "covariates only", repeat_id = r, fold = k,
            metric = names(m), value = unname(m), n_test = nrow(test))
        }
        for (nm in names(specs)) {
          n_valid <- nrow(test)
          res <- tryCatch({
            ref <- if (is.null(reference_data)) fit_reference(train, specs[[nm]]) else fit_reference(reference_data, specs[[nm]])
            s_train <- score_load(train, ref)$scores$score; s_test <- score_load(test, ref)$scores$score
            ok_train <- !is.na(s_train); ok_test <- !is.na(s_test)
            n_valid <- sum(ok_test)
            fit_predict(train[ok_train, , drop = FALSE], test[ok_test, , drop = FALSE], s_train[ok_train], s_test[ok_test])
          }, error = function(e) {
            failures[[length(failures) + 1L]] <<- data.frame(spec = nm, repeat_id = r, fold = k, error = conditionMessage(e),
              stringsAsFactors = FALSE)
            NULL
          })
          if (!is.null(res)) out[[length(out) + 1L]] <- data.frame(spec = nm, repeat_id = r, fold = k, metric = names(res),
            value = unname(res), n_test = n_valid)
        }
      }
    }
    if (length(out)) do.call(rbind, out) else data.frame(spec = character(), repeat_id = integer(), fold = integer(),
      metric = character(), value = numeric(), n_test = integer())
  })
  failures_df <- if (length(failures)) do.call(rbind, failures) else data.frame(spec = character(), repeat_id = integer(),
    fold = integer(), error = character())
  if (nrow(failures_df)) {
    failed_specs <- unique(failures_df$spec)
    warning(sprintf("Cross-validation failed for %s (%d fold(s)): %s", paste(failed_specs, collapse = ", "),
      nrow(failures_df), failures_df$error[1]), call. = FALSE)
  }
  if (!nrow(rows)) stop("All specifications failed during cross-validation.", call. = FALSE)
  summary <- stats::aggregate(value ~ spec + metric, data = rows, FUN = function(v) c(mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE)))
  summary <- data.frame(spec = summary$spec, metric = summary$metric, mean = summary$value[, "mean"], sd = summary$value[, "sd"])
  main <- switch(type, binary = "auc", survival = "c_index", continuous = "r_squared")
  if (length(covariates)) {
    base <- summary$mean[summary$spec == "covariates only" & summary$metric == main]
    summary$increment_over_covariates <- ifelse(summary$metric == main & summary$spec != "covariates only", summary$mean - base, NA_real_)
  }
  lower_better <- summary$metric %in% c("brier", "rmse")
  summary <- summary[order(summary$metric, ifelse(lower_better, summary$mean, -summary$mean)), ]; rownames(summary) <- NULL
  shown <- summary[summary$metric == main, ]
  shown$spec <- factor(shown$spec, levels = shown$spec[order(shown$mean)])
  label <- switch(main, auc = "Cross-validated AUC", c_index = "Cross-validated Harrell's C", r_squared = "Cross-validated R-squared")
  plot <- .finish_plot(
    ggplot2::ggplot(shown, ggplot2::aes(x = .data$mean, y = .data$spec, xmin = .data$mean - .data$sd, xmax = .data$mean + .data$sd)) +
      ggplot2::geom_pointrange(colour = .al_palette$groups[["ink"]]),
    title = "Out-of-sample performance of scoring algorithms", x = label, y = NULL, width = 100,
    caption = sprintf("%d-fold cross-validation x %d repeat(s); bars: +/- 1 SD across folds. References %s.",
      folds, repeats, if (is.null(reference_data)) "refitted in every training fold" else "fixed (external)"))
  structure(list(folds = rows, summary = summary, failures = failures_df, type = type, outcome = outcome,
    plots = list(performance = plot)), class = "al_validation")
}

#' @export
print.al_validation <- function(x, ...) {
  cat(sprintf("Cross-validated validation against %s (%s outcome)\n\n", x$outcome, x$type))
  print(x$summary, row.names = FALSE, digits = 3, ...)
  invisible(x)
}
