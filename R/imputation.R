.imputation_methods <- c("rf", "missforest", "pmm", "cart", "norm")

.prepare_imputation <- function(data, spec, predictors) {
  validate_data(data, spec)
  names_all <- .biomarker_field(spec, "name")
  if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors))
    stop("predictors must be a character vector of unique column names.", call. = FALSE)
  predictors <- unique(c(predictors, spec$strata, spec$age))
  if (any(predictors %in% names_all)) stop("Biomarkers are imputed, not listed as predictors.", call. = FALSE)
  if (!all(predictors %in% names(data))) stop("Missing predictor columns: ", paste(setdiff(predictors, names(data)), collapse = ", "), call. = FALSE)
  for (p in predictors) {
    x <- data[[p]]
    if (!is.atomic(x) || !is.null(dim(x)) || anyNA(x)) stop("Predictor ", p, " must be a complete vector.", call. = FALSE)
  }
  for (b in names_all) {
    if (all(is.na(data[[b]]))) stop("Biomarker '", b, "' has no observed values and cannot be imputed.", call. = FALSE)
  }
  transforms <- .biomarker_field(spec, "transform")
  work <- data.frame(lapply(seq_along(names_all), function(j) .transform(data[[names_all[j]]], transforms[j])))
  names(work) <- sprintf("b%02d", seq_along(names_all))
  for (k in seq_along(predictors)) {
    x <- data[[predictors[k]]]
    work[[sprintf("p%02d", k)]] <- if (is.character(x) || is.logical(x)) factor(x) else x
  }
  list(work = work, names = names_all, transforms = transforms, predictors = predictors,
    mask = is.na(as.matrix(data[names_all])))
}

# Native missForest (Stekhoven and Buhlmann 2012): iterate random forests until
# the imputed values stop changing; return the last imputation before the
# change increased.
.missforest <- function(work, columns, maxit, ntree) {
  missing <- is.na(work[columns])
  x <- work
  for (v in columns) x[[v]][is.na(x[[v]])] <- mean(x[[v]], na.rm = TRUE)
  order <- columns[order(colSums(missing))]
  previous <- x; difference <- Inf; oob <- stats::setNames(rep(NA_real_, length(columns)), columns)
  iterations <- 0L
  repeat {
    iterations <- iterations + 1L
    current <- previous; current_oob <- oob
    for (v in order) {
      observed <- !missing[, v]
      if (all(observed)) next
      fit <- randomForest::randomForest(x = current[observed, setdiff(names(current), v), drop = FALSE],
        y = current[[v]][observed], ntree = ntree)
      current[[v]][!observed] <- stats::predict(fit, current[!observed, setdiff(names(current), v), drop = FALSE])
      current_oob[v] <- sqrt(fit$mse[ntree] / stats::var(current[[v]][observed]))
    }
    new_difference <- sum((as.matrix(current[columns]) - as.matrix(previous[columns]))^2) /
      sum(as.matrix(current[columns])^2)
    if (new_difference >= difference) { iterations <- iterations - 1L; break }
    previous <- current; oob <- current_oob; difference <- new_difference
    if (iterations >= maxit) break
  }
  list(data = previous, oob = oob, iterations = iterations)
}

#' Impute missing biomarker values
#'
#' Explicit, auditable imputation of missing biomarker measurements before
#' scoring. Nothing is imputed unless this function is called; the original
#' data are never modified and every imputed cell is recorded.
#'
#' `"rf"` draws multiple imputations from random forests within chained
#' equations (mice's `rf` method with the randomForest package; Doove et al.
#' 2014). `"missforest"` implements the iterative random-forest algorithm of
#' Stekhoven and Buhlmann (2012) and returns a single imputation with
#' out-of-bag error estimates; single imputation understates uncertainty.
#' `"pmm"` (predictive mean matching), `"cart"` and `"norm"` (Bayesian linear
#' regression) use mice. Biomarkers with a log transform are imputed on the
#' log scale and back-transformed, so imputed values stay positive.
#' Protocol strata and age columns are always used as predictors.
#'
#' @param data Data frame with protocol biomarkers (missing values allowed).
#' @inheritParams validate_data
#' @param method Imputation method; see Details.
#' @param m Number of imputations (1 for `"missforest"`).
#' @param predictors Additional complete columns used as predictors, never imputed.
#' @param maxit Maximum number of iterations.
#' @param ntree Number of trees per forest; defaults to 10 for `"rf"` (mice's
#'   default) and 100 for `"missforest"`.
#' @inheritParams simulate_biomarkers
#' @param delta Optional named numeric shifts for a delta-adjustment
#'   sensitivity analysis of data missing not at random: added to imputed
#'   values only, on the imputation scale (a log-scale shift multiplies
#'   log-transformed biomarkers by `exp(delta)`). See [tipping_point()].
#' @return An `al_imputation` object with `datasets` (m completed copies of
#'   `data`), `mask` (TRUE for imputed cells), `diagnostics` (missingness per
#'   biomarker, chain means for mice methods, out-of-bag normalized RMSE for
#'   missforest and mice logged events), software versions and plots
#'   (`distributions`; `convergence` for mice methods). Score it with
#'   [score_imputed()] and fit pooled references with [fit_reference()].
#' @references van Buuren S, Groothuis-Oudshoorn K (2011). mice: Multivariate
#'   imputation by chained equations in R. Journal of Statistical Software 45(3).
#'   Stekhoven DJ, Buhlmann P (2012). MissForest: non-parametric missing value
#'   imputation for mixed-type data. Bioinformatics 28, 112-118.
#'   Doove LL, van Buuren S, Dusseldorp E (2014). Recursive partitioning for
#'   missing data imputation in the presence of interaction effects.
#'   Computational Statistics and Data Analysis 72, 92-104.
#' @family imputation
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
impute_biomarkers <- function(data, spec, method = c("rf", "missforest", "pmm", "cart", "norm"), m = 5L,
                              predictors = character(), maxit = 10L, ntree = NULL, seed = 1L, delta = numeric()) {
  method <- match.arg(method)
  if (!is.numeric(delta) || any(!is.finite(delta)) || (length(delta) && (!.valid_names(names(delta)) ||
      !all(names(delta) %in% .biomarker_field(spec, "name")))))
    stop("delta must be a named numeric vector of biomarker shifts.", call. = FALSE)
  .check_integer(m, "m", 1L, 100L); .check_integer(maxit, "maxit", 1L, 100L)
  if (method == "missforest" && m != 1L) stop("missforest produces a single imputation; use m = 1 or method = 'rf'.", call. = FALSE)
  if (is.null(ntree)) ntree <- if (method == "missforest") 100L else 10L
  .check_integer(ntree, "ntree", 1L)
  if (method %in% c("rf", "missforest")) .require_suggested("randomForest", "impute with random forests")
  if (method != "missforest") .require_suggested("mice", "impute with chained equations")
  prep <- .prepare_imputation(data, spec, predictors)
  work <- prep$work; biomarker_columns <- sprintf("b%02d", seq_along(prep$names))
  columns <- biomarker_columns[colSums(prep$mask) > 0]
  diagnostics <- list(missing = data.frame(biomarker = prep$names, n_missing = colSums(prep$mask),
    fraction = colMeans(prep$mask), row.names = NULL), all_missing_rows = sum(rowSums(prep$mask) == ncol(prep$mask)))
  completed <- list()
  if (!length(columns)) {
    completed <- rep(list(work), m)
  } else if (method == "missforest") {
    fit <- .with_seed(seed, function() .missforest(work, columns, maxit, ntree))
    completed <- list(fit$data)
    diagnostics$oob <- data.frame(biomarker = prep$names[match(columns, biomarker_columns)], nrmse = unname(fit$oob[columns]))
    diagnostics$iterations <- fit$iterations
  } else {
    methods <- stats::setNames(ifelse(names(work) %in% columns, method, ""), names(work))
    extra <- if (method == "rf") list(ntree = ntree, rfPackage = "randomForest") else list()
    imp <- .with_seed(seed, function() do.call(mice::mice, c(list(data = work, m = m, maxit = maxit,
      method = methods, printFlag = FALSE), extra)))
    completed <- lapply(seq_len(m), function(k) mice::complete(imp, action = k))
    chain <- imp$chainMean
    diagnostics$chains <- do.call(rbind, lapply(columns, function(v) {
      values <- chain[v, , , drop = FALSE]
      data.frame(biomarker = prep$names[match(v, biomarker_columns)],
        iteration = rep(seq_len(dim(chain)[2]), times = dim(chain)[3]),
        chain = rep(seq_len(dim(chain)[3]), each = dim(chain)[2]), mean = as.vector(values))
    }))
    diagnostics$logged_events <- imp$loggedEvents
  }
  datasets <- lapply(completed, function(w) {
    out <- data
    for (j in seq_along(prep$names)) {
      missing <- prep$mask[, j]
      shift <- if (prep$names[j] %in% names(delta)) delta[[prep$names[j]]] else 0
      if (any(missing)) out[[prep$names[j]]][missing] <- .untransform(w[[biomarker_columns[j]]][missing] + shift, prep$transforms[j])
    }
    out
  })
  packages <- c("allostaticR", if (method != "missforest") "mice", if (method %in% c("rf", "missforest")) "randomForest")
  out <- structure(list(datasets = datasets, mask = prep$mask, method = method, m = m, maxit = maxit,
    ntree = if (method %in% c("rf", "missforest")) ntree else NA_integer_, seed = seed,
    predictors = prep$predictors, imputed = prep$names[colSums(prep$mask) > 0], spec = spec, delta = delta,
    diagnostics = diagnostics, packages = data.frame(package = packages,
      version = vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), row.names = NULL),
    plots = list()), class = "al_imputation")
  out$plots <- .imputation_plots(out, data)
  out
}

.imputation_plots <- function(x, data) {
  if (!length(x$imputed)) return(list())
  info <- .biomarker_info(x$spec); colors <- .system_colors(.protocol_systems(x$spec))
  points <- do.call(rbind, lapply(x$imputed, function(b) {
    j <- match(b, info$name); observed <- data[[b]][!x$mask[, j]]
    if (length(observed) > 1000L) observed <- .with_seed(1L, function() sample(observed, 1000L))
    imputed <- do.call(rbind, lapply(seq_along(x$datasets), function(k) {
      vals <- x$datasets[[k]][[b]][x$mask[, j]]
      if (!length(vals)) return(NULL)
      data.frame(source = paste("Imp.", k), value = vals)
    }))
    obs_df <- if (length(observed)) data.frame(source = "Observed", value = observed) else data.frame(source = character(), value = numeric())
    d <- rbind(obs_df, imputed)
    d$biomarker <- info$label[j]; d$system <- info$system[j]
    d$kind <- ifelse(d$source == "Observed", "Observed", "Imputed")
    d
  }))
  points$source <- factor(points$source, c("Observed", paste("Imp.", seq_along(x$datasets))))
  points$biomarker <- factor(points$biomarker, info$label[info$name %in% x$imputed])
  plots <- list(distributions = .finish_plot(
    ggplot2::ggplot(points, ggplot2::aes(x = .data$source, y = .data$value)) +
      ggplot2::geom_point(data = points[points$kind == "Observed", ], colour = .al_palette$light_grey, alpha = 0.35, size = 0.8,
        position = ggplot2::position_jitter(width = 0.25, height = 0, seed = 1)) +
      ggplot2::geom_point(data = points[points$kind == "Imputed", ], ggplot2::aes(colour = .data$system), alpha = 0.8,
        size = 1.1, position = ggplot2::position_jitter(width = 0.25, height = 0, seed = 2)) +
      ggplot2::scale_colour_manual(values = colors, breaks = names(colors)) +
      ggplot2::facet_wrap(~biomarker, scales = "free_y"),
    title = "Observed and imputed values", x = NULL, y = "Value (declared unit)", colour = "Imputed",
    subtitle = sprintf("Method: %s | imputations: %d | gray: observed values", x$method, x$m),
    caption = "Imputed values should be plausible relative to observed ones; MAR is assumed, not verified.") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  if (!is.null(x$diagnostics$chains) && x$maxit > 1L) {
    chains <- x$diagnostics$chains
    plots$convergence <- .finish_plot(
      ggplot2::ggplot(chains, ggplot2::aes(x = .data$iteration, y = .data$mean, group = .data$chain)) +
        ggplot2::geom_line(colour = .al_palette$secondary, alpha = 0.7) +
        ggplot2::facet_wrap(~biomarker, scales = "free_y"),
      title = "Convergence of imputation chains", x = "Iteration",
      y = "Mean of imputed values (imputation scale)",
      caption = "Chains should intermingle without trends; log-transformed biomarkers are shown on the log scale.")
  }
  plots
}

#' @export
print.al_imputation <- function(x, ...) {
  cat(sprintf("Biomarker imputation: %s | imputations: %d | imputed cells: %d\n", x$method, x$m, sum(x$mask)))
  cat("Predictors:", if (length(x$predictors)) paste(x$predictors, collapse = ", ") else "biomarkers only", "\n\n")
  print(x$diagnostics$missing, row.names = FALSE, digits = 3, ...)
  if (!is.null(x$diagnostics$oob)) { cat("\nOut-of-bag normalized RMSE:\n"); print(x$diagnostics$oob, row.names = FALSE, digits = 3) }
  if (NROW(x$diagnostics$logged_events)) cat("\nmice logged", NROW(x$diagnostics$logged_events),
    "events (for example predictors removed for collinearity); see diagnostics$logged_events.\n")
  invisible(x)
}

.fit_imputed_reference <- function(imputation, spec, label) {
  refs <- lapply(imputation$datasets, fit_reference, spec = spec, label = label)
  ref <- refs[[1]]
  for (column in intersect(c("lower", "upper", "center", "scale"), names(ref$cutoffs)))
    ref$cutoffs[[column]] <- rowMeans(vapply(refs, function(r) r$cutoffs[[column]], numeric(nrow(ref$cutoffs))))
  if (!is.null(ref$models)) for (k in seq_along(ref$models)) for (p in names(ref$models[[k]]$coef))
    ref$models[[k]]$coef[[p]] <- rowMeans(vapply(refs, function(r) r$models[[k]]$coef[[p]], numeric(length(ref$models[[k]]$coef[[p]]))))
  ref$imputation <- list(method = imputation$method, m = imputation$m, imputed_cells = sum(imputation$mask))
  .assert_reference(ref)
  ref
}

# Rubin's rules with the Barnard-Rubin small-sample degrees of freedom.
.rubin <- function(estimates, variances, df_complete, conf = 0.95) {
  m <- length(estimates); qbar <- mean(estimates); ubar <- mean(variances)
  b <- if (m > 1L) stats::var(estimates) else 0
  total <- ubar + (1 + 1 / m) * b
  lambda <- if (total > 0) (1 + 1 / m) * b / total else 0
  df_observed <- (df_complete + 1) / (df_complete + 3) * df_complete * (1 - lambda)
  df <- if (b > 0) 1 / (lambda^2 / (m - 1) + 1 / df_observed) else df_complete
  riv <- if (ubar > 0) (1 + 1 / m) * b / ubar else NA_real_
  half <- stats::qt(1 - (1 - conf) / 2, df) * sqrt(total)
  data.frame(estimate = qbar, std_error = sqrt(total), lower = qbar - half, upper = qbar + half, df = df,
    within = ubar, between = b, lambda = lambda, fmi = if (is.finite(riv)) (riv + 2 / (df + 3)) / (riv + 1) else NA_real_, m = m)
}

#' Score imputed datasets and pool the results
#'
#' Scores every completed dataset of an imputation against the same fixed
#' reference and pools cohort means (overall, by group and the difference
#' between two groups) with Rubin's rules and Barnard-Rubin degrees of
#' freedom.
#'
#' @param imputation An `al_imputation` from impute_biomarkers().
#' @param reference Fixed reference from fit_reference().
#' @param id,group Optional identifier and group columns.
#' @param levels Optional two groups for the pooled difference (`levels[2]`
#'   minus `levels[1]`); defaults to the two groups when exactly two exist.
#' @param conf Confidence level of pooled intervals.
#' @return An `al_imputed_result` with `results` (one al_result per
#'   imputation, usable with every other function), `scores` (per
#'   participant: mean, minimum and maximum score across imputations and the
#'   number of imputed biomarkers), `components` (share of imputations in
#'   which each biomarker is at risk), `pooled` (Rubin-pooled means and the
#'   group difference) and a plot of pooled group means. With a single
#'   imputation the between-imputation variance is zero and intervals are too narrow.
#' @references Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys.
#'   Wiley. Barnard J, Rubin DB (1999). Small-sample degrees of freedom with
#'   multiple imputation. Biometrika 86, 948-955.
#' @family imputation
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
score_imputed <- function(imputation, reference, id = NULL, group = NULL, levels = NULL, conf = 0.95) {
  if (!inherits(imputation, "al_imputation")) stop("imputation must come from impute_biomarkers().", call. = FALSE)
  .assert_reference(reference); .check_open_unit(conf, "conf")
  if (!all(.biomarker_field(reference$spec, "name") %in% colnames(imputation$mask)))
    stop("reference biomarkers must belong to the imputed panel.", call. = FALSE)
  results <- lapply(imputation$datasets, score_load, reference = reference, id = id, group = group)
  first <- results[[1]]$scores
  score_mat <- vapply(results, function(r) r$scores$score, numeric(nrow(first)))
  score_mat <- matrix(score_mat, nrow = nrow(first))
  imputed_names <- colnames(imputation$mask)
  scores <- data.frame(row = first$row, id = first$id, score_mean = rowMeans(score_mat),
    score_min = apply(score_mat, 1, min), score_max = apply(score_mat, 1, max),
    n_imputed = rowSums(imputation$mask[, .biomarker_field(reference$spec, "name"), drop = FALSE]),
    status = ifelse(rowSums(is.na(score_mat)) == 0, "ok", ifelse(rowSums(is.na(score_mat)) == ncol(score_mat),
      first$status, "invalid_in_some_imputations")))
  if (!is.null(group)) scores$group <- first$group
  contribution <- vapply(results, function(r) r$components$contribution, numeric(nrow(results[[1]]$components)))
  contribution <- matrix(contribution, nrow = nrow(results[[1]]$components))
  cmp <- results[[1]]$components
  components <- data.frame(row = cmp$row, id = cmp$id, biomarker = cmp$biomarker,
    imputed = as.vector(imputation$mask[, match(unique(cmp$biomarker), imputed_names)]),
    contribution_mean = rowMeans(contribution))
  if (.is_count_method(reference$spec)) components$risk_share <- components$contribution_mean
  pool_mean <- function(rows, label) {
    q <- colMeans(score_mat[rows, , drop = FALSE]); u <- apply(score_mat[rows, , drop = FALSE], 2, stats::var) / length(rows)
    cbind(target = label, n = length(rows), .rubin(q, u, length(rows) - 1, conf))
  }
  valid <- rowSums(is.na(score_mat)) == 0
  if (sum(valid) < 2L) stop("At least two participants must be scored in every imputation.", call. = FALSE)
  pooled <- list(pool_mean(which(valid), "all"))
  if (!is.null(group)) {
    g <- as.character(first$group)
    for (lv in sort(unique(stats::na.omit(g)))) if (sum(valid & g %in% lv) > 1L) pooled[[length(pooled) + 1L]] <- pool_mean(which(valid & g %in% lv), lv)
    available <- sort(unique(stats::na.omit(g[valid])))
    if (is.null(levels) && length(available) == 2L) levels <- available
    if (!is.null(levels)) {
      .check_group_levels(levels, available, 2L)
      a <- which(valid & g %in% levels[1]); b <- which(valid & g %in% levels[2])
      q <- colMeans(score_mat[b, , drop = FALSE]) - colMeans(score_mat[a, , drop = FALSE])
      u <- apply(score_mat[a, , drop = FALSE], 2, stats::var) / length(a) + apply(score_mat[b, , drop = FALSE], 2, stats::var) / length(b)
      pooled[[length(pooled) + 1L]] <- cbind(target = paste(levels[2], "minus", levels[1]), n = length(a) + length(b),
        .rubin(q, u, length(a) + length(b) - 2, conf))
    }
  }
  pooled <- do.call(rbind, pooled); rownames(pooled) <- NULL
  means <- pooled[!grepl(" minus ", pooled$target), ]
  plot <- .finish_plot(
    ggplot2::ggplot(means, ggplot2::aes(x = .data$estimate, y = factor(.data$target, rev(.data$target)),
      xmin = .data$lower, xmax = .data$upper)) +
      ggplot2::geom_pointrange(colour = .al_palette$groups[["ink"]]),
    title = "Pooled mean scores after multiple imputation", x = paste("Mean", tolower(.scale_label(results[[1]]$audit$scale))),
    y = NULL, width = 100,
    caption = sprintf("Rubin's rules across %d imputations (%s); %g%% intervals with Barnard-Rubin df.",
      imputation$m, imputation$method, 100 * conf))
  structure(list(results = results, scores = scores, components = components, pooled = pooled,
    imputation = imputation[c("method", "m", "seed", "predictors", "imputed")], plots = list(pooled_means = plot)),
    class = "al_imputed_result")
}

#' @export
print.al_imputed_result <- function(x, ...) {
  cat(sprintf("Scores pooled over %d imputations (%s)\n\n", x$imputation$m, x$imputation$method))
  print(x$pooled[c("target", "n", "estimate", "lower", "upper", "df", "fmi")], row.names = FALSE, digits = 3, ...)
  invisible(x)
}

#' Compare imputation methods by masking observed values
#'
#' Masks a random share of observed biomarker values (missing completely at
#' random), imputes them with each method and compares the imputations with
#' the hidden true values: normalized RMSE, bias and, for count protocols,
#' agreement of the resulting risk flags with the true flags. A median
#' baseline is always included. Results describe accuracy under the masking
#' mechanism only; real missingness may be informative.
#'
#' @param data Data frame with protocol biomarkers.
#' @inheritParams validate_data
#' @param reference Fixed reference used to compare risk flags.
#' @param methods Imputation methods to compare.
#' @param prop Share of observed values masked per biomarker.
#' @param m,maxit,ntree,predictors Passed to impute_biomarkers().
#' @param seed Local random seed for masking and imputation.
#' @return An `al_imputation_check` with `accuracy` (per method and biomarker),
#'   `summary` (per method) and editable plots `nrmse` and `flags`.
#' @family imputation
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
check_imputation <- function(data, spec, reference, methods = c("rf", "pmm", "cart", "norm"), prop = 0.1,
                             m = 5L, maxit = 5L, ntree = NULL, predictors = character(), seed = 1L) {
  methods <- match.arg(methods, .imputation_methods, several.ok = TRUE)
  .assert_reference(reference); .check_open_unit(prop, "prop")
  if (!identical(reference$spec$biomarkers, spec$biomarkers)) stop("reference must belong to the same biomarker panel.", call. = FALSE)
  prep <- .prepare_imputation(data, spec, predictors)
  names_all <- prep$names
  masked_cells <- .with_seed(seed, function() lapply(seq_along(names_all), function(j) {
    observed <- which(!prep$mask[, j]); k <- max(1L, round(prop * length(observed)))
    if (length(observed) < 3L) stop("Too few observed values to mask for ", names_all[j], ".", call. = FALSE)
    observed[sample.int(length(observed), k)]
  }))
  masked <- data
  for (j in seq_along(names_all)) masked[[names_all[j]]][masked_cells[[j]]] <- NA
  params <- .reference_parameters(reference, data)
  count <- .is_count_method(spec)
  rows <- list()
  for (method in c("median", methods)) {
    if (method == "median") {
      baseline <- masked
      for (b in names_all) baseline[[b]][is.na(baseline[[b]])] <- stats::median(baseline[[b]], na.rm = TRUE)
      datasets <- list(baseline)
    } else datasets <- impute_biomarkers(masked, spec, method = method, m = if (method == "missforest") 1L else m,
      predictors = predictors, maxit = maxit, ntree = ntree, seed = seed)$datasets
    for (j in seq_along(names_all)) {
      cells <- masked_cells[[j]]; truth <- data[[names_all[j]]][cells]
      imputed <- matrix(vapply(datasets, function(d) d[[names_all[j]]][cells], numeric(length(cells))), nrow = length(cells))
      point <- rowMeans(imputed)
      b <- spec$biomarkers[[j]]; p <- lapply(params[[j]][c("lower", "upper")], `[`, cells)
      agreement <- kappa <- NA_real_
      if (count) {
        true_flag <- .risk_flag(truth, p, b$direction, spec$inclusive)
        flags <- apply(imputed, 2, .risk_flag, p = p, direction = b$direction, inclusive = spec$inclusive)
        flags <- matrix(flags, nrow = length(cells))
        agreement <- mean(flags == true_flag, na.rm = TRUE)
        kappa <- .kappa(rowMeans(flags) >= 0.5, true_flag)
      }
      rows[[length(rows) + 1L]] <- data.frame(method = method, biomarker = names_all[j], n_masked = length(cells),
        nrmse = sqrt(mean((point - truth)^2)) / stats::sd(data[[names_all[j]]], na.rm = TRUE),
        bias = mean(point - truth), flag_agreement = agreement, kappa = kappa)
    }
  }
  accuracy <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(accuracy, factor(accuracy$method, unique(accuracy$method))), function(d)
    data.frame(method = d$method[1], mean_nrmse = mean(d$nrmse), flag_agreement = stats::weighted.mean(d$flag_agreement, d$n_masked),
      mean_kappa = .mean_or_na(d$kappa))))
  rownames(summary) <- NULL
  info <- .biomarker_info(spec)
  shown <- accuracy
  shown$label <- factor(info$label[match(shown$biomarker, info$name)], levels = info$label)
  shown$method <- factor(shown$method, rev(unique(shown$method)))
  baseline <- shown[shown$method == "median", ]
  panel <- function(column, title, x_label, percent) {
    .finish_plot(
      ggplot2::ggplot(shown[shown$method != "median", ], ggplot2::aes(x = .data[[column]], y = .data$method)) +
        ggplot2::geom_vline(data = baseline, ggplot2::aes(xintercept = .data[[column]]), colour = .al_palette$muted,
          linetype = "22", linewidth = 0.5) +
        ggplot2::geom_point(size = 2.4, colour = .al_palette$secondary) +
        (if (percent) ggplot2::scale_x_continuous(labels = .percent) else NULL) +
        ggplot2::facet_wrap(~label, ncol = 3),
      title = title, x = x_label, y = NULL, width = 100,
      caption = sprintf("%s of observed values masked completely at random per biomarker. Dashed line: median imputation baseline.",
        .percent(prop)))
  }
  plots <- list(nrmse = panel("nrmse", "Imputation accuracy on masked values", "Normalized RMSE (lower is better)", FALSE))
  if (count) plots$flags <- panel("flag_agreement", "Risk flags recovered after imputation", "Agreement with the true risk flag", TRUE)
  accuracy$method <- factor(accuracy$method, unique(accuracy$method))
  structure(list(accuracy = accuracy, summary = summary, prop = prop, plots = plots), class = "al_imputation_check")
}

#' @export
print.al_imputation_check <- function(x, ...) {
  cat(sprintf("Imputation check: %s of observed values masked per biomarker\n\n", .percent(x$prop)))
  print(x$summary, row.names = FALSE, digits = 3, ...)
  invisible(x)
}

#' Tipping-point analysis for data missing not at random
#'
#' Repeats multiple imputation while shifting the imputed values of one
#' biomarker by a range of deltas (delta adjustment), then pools the group
#' difference in mean score at each delta. The tipping point is the smallest
#' shift at which the conclusion changes (the interval first includes zero).
#' Plausible deltas must be justified substantively.
#'
#' @param data,spec,reference Data, protocol and fixed reference.
#' @param biomarker Biomarker whose imputed values are shifted.
#' @param deltas Shifts on the imputation scale (log scale for log-transformed
#'   biomarkers); include 0 for the missing-at-random analysis.
#' @param group Group column; `levels` chooses the comparison.
#' @param levels Optional two groups (`levels[2]` minus `levels[1]`).
#' @param method,m,maxit,predictors,seed Passed to impute_biomarkers(); the
#'   same seed is used at every delta.
#' @return An `al_tipping_point` with `results` (pooled difference per delta),
#'   `tipping_delta` (first delta, in order of distance from 0, whose interval
#'   includes 0, or NA) and an editable plot.
#' @family imputation
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
tipping_point <- function(data, spec, reference, biomarker, deltas, group, levels = NULL, method = "pmm",
                          m = 5L, maxit = 5L, predictors = character(), seed = 1L) {
  .scalar_string(biomarker, "biomarker"); .scalar_string(group, "group")
  if (!is.numeric(deltas) || !length(deltas) || any(!is.finite(deltas))) stop("deltas must be finite numbers.", call. = FALSE)
  if (!sum(is.na(data[[biomarker]]))) stop("The biomarker has no missing values to shift.", call. = FALSE)
  deltas <- sort(unique(deltas))
  results <- do.call(rbind, lapply(deltas, function(d) {
    imp <- impute_biomarkers(data, spec, method = method, m = m, maxit = maxit, predictors = predictors, seed = seed,
      delta = stats::setNames(d, biomarker))
    pooled <- score_imputed(imp, reference, group = group, levels = levels)$pooled
    row <- pooled[grepl(" minus ", pooled$target), ]
    if (!nrow(row)) stop("Supply levels for the group comparison.", call. = FALSE)
    data.frame(delta = d, target = row$target, estimate = row$estimate, lower = row$lower, upper = row$upper, fmi = row$fmi)
  }))
  includes <- results$lower <= 0 & results$upper >= 0
  base_excludes <- !includes[which.min(abs(results$delta))]
  order_away <- order(abs(results$delta))
  tipping <- if (base_excludes && any(includes)) results$delta[order_away][includes[order_away]][1] else NA_real_
  unit <- if (spec$biomarkers[[match(biomarker, .biomarker_field(spec, "name"))]]$transform == "log")
    "log-scale shift of imputed values" else "shift of imputed values (declared unit)"
  plot <- .finish_plot(
    ggplot2::ggplot(results, ggplot2::aes(x = .data$delta, y = .data$estimate)) +
      ggplot2::geom_hline(yintercept = 0, colour = .al_palette$secondary, linewidth = 0.4) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), fill = .al_palette$light_blue) +
      ggplot2::geom_line(colour = .al_palette$categorical[[1]], linewidth = 0.9) + ggplot2::geom_point(colour = .al_palette$categorical[[1]], size = 2) +
      (if (is.finite(tipping)) ggplot2::geom_vline(xintercept = tipping, linetype = "22", colour = .al_palette$diverging[["high"]]) else NULL),
    title = paste("Tipping-point analysis:", biomarker),
    subtitle = if (is.finite(tipping)) sprintf("Conclusion changes at delta = %g", tipping) else "Conclusion unchanged across the deltas examined",
    x = paste("Delta:", unit), y = paste("Pooled difference:", results$target[1]), width = 100,
    caption = sprintf("Delta-adjusted multiple imputation (%s, m = %d, same seed at every delta). Delta 0 assumes missing at random.", method, m))
  structure(list(results = results, tipping_delta = tipping, biomarker = biomarker, plots = list(tipping = plot)),
    class = "al_tipping_point")
}

#' @export
print.al_tipping_point <- function(x, ...) {
  cat("Tipping-point analysis for", x$biomarker, "| tipping delta:", if (is.finite(x$tipping_delta)) x$tipping_delta else "none", "\n\n")
  print(x$results, row.names = FALSE, digits = 3, ...)
  invisible(x)
}
