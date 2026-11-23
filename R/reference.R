.transform <- function(x, transform) if (identical(transform, "log")) log(x) else x
.untransform <- function(x, transform) if (identical(transform, "log")) exp(x) else x

.stratum_key <- function(data, strata) {
  if (!length(strata)) return(rep("all", nrow(data)))
  do.call(paste, c(lapply(strata, function(s) as.character(data[[s]])), sep = " | "))
}

# Measured values, with fixed treatment effects added under medication_rule = "add".
.adjusted_values <- function(data, spec, b) {
  x <- data[[b$name]]
  if (spec$medication_rule != "add" || !b$name %in% names(spec$medications)) return(x)
  column <- spec$medications[[b$name]]
  if (length(.invalid_medication_columns(data, column)))
    stop("Medication columns must be complete logical vectors: ", column, call. = FALSE)
  x + spec$medication_effects[[b$name]] * data[[column]]
}

.reference_weights <- function(data, spec) {
  if (!length(spec$weights)) return(rep(1, nrow(data)))
  if (!spec$weights %in% names(data)) stop("Missing weights column: ", spec$weights, call. = FALSE)
  w <- data[[spec$weights]]
  if (!is.numeric(w) || any(!is.finite(w)) || any(w <= 0)) stop("Weights must be positive and finite.", call. = FALSE)
  w
}

# Weighted quantile by linear interpolation of mid-point cumulative weights
# (equals quantile type 5 with equal weights).
.weighted_quantile <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  position <- (cumsum(w) - w / 2) / sum(w)
  stats::approx(position, x, xout = p, rule = 2, ties = "ordered")$y
}

.in_stratum <- function(level, spec) if (length(spec$strata)) paste0(" in stratum '", level, "'") else ""

# Tail probabilities defining quantile thresholds for one biomarker.
.tail_probs <- function(spec, direction) {
  rq <- spec$risk_quantile
  if (direction == "high") return(c(upper = rq))
  if (direction == "low") return(c(lower = 1 - rq))
  if (spec$bilateral == "split") c(lower = (1 - rq) / 2, upper = 1 - (1 - rq) / 2) else c(lower = 1 - rq, upper = rq)
}

.spline_design <- function(age, model) {
  cbind(1, splines::ns(age, knots = model$knots, Boundary.knots = model$boundary))
}

.cutoff_row <- function(b, level, spec, lower = NA_real_, upper = NA_real_, center = NA_real_,
                        scale = NA_real_, n = NA_integer_, age_range = NULL) {
  d <- data.frame(biomarker = b$name, label = b$label, system = b$system, direction = b$direction,
    unit = b$unit)
  if (length(spec$strata)) d$stratum <- level
  d$lower <- lower; d$upper <- upper
  if (spec$method == "z_mean") { d$center <- center; d$scale <- scale }
  d$n_reference <- as.integer(n)
  if (length(spec$age)) { d$age_min <- age_range[1]; d$age_max <- age_range[2] }
  d
}

.new_reference <- function(spec, cutoffs, label, n_reference, levels, strata_values = NULL, models = NULL) {
  rownames(cutoffs) <- NULL
  structure(list(spec = spec, cutoffs = cutoffs, label = label, n_reference = n_reference,
    levels = levels, strata_values = strata_values, models = models), class = "al_reference")
}

.strata_values <- function(data, spec, levels) {
  if (!length(spec$strata)) return(NULL)
  key <- .stratum_key(data, spec$strata)
  d <- data[match(levels, key), spec$strata, drop = FALSE]
  rownames(d) <- NULL
  d
}

#' Fit reusable biomarker thresholds
#'
#' Estimates the fixed reference parameters of a protocol: quantile thresholds
#' (`quantile_count`), supplied cutoffs (`cutoff_count`) or reference means and
#' SDs (`z_mean`). Protocol strata are fitted separately; an age column yields
#' smooth age-specific thresholds (quantile regression) or age-specific means
#' and SDs (Altman's absolute-residual method), like growth-chart centiles.
#'
#' @param reference Reference data frame. Optional for supplied cutoffs. An
#'   `al_imputation` object from [impute_biomarkers()] is also accepted: the
#'   reference is fitted in each completed dataset and the parameters are averaged.
#' @inheritParams validate_data
#' @param label Name identifying the reference population or cutoff protocol.
#' @return An `al_reference` object with fixed thresholds and reference sample sizes.
#'   Age models need at least `10 * (age_df + 1)` observed values per biomarker
#'   and stratum, and thresholds that cross are rejected.
#' @references Altman DG (1993). Construction of age-related reference centiles
#'   using absolute residuals. Statistics in Medicine 12, 917-924.
#' @family reference calibration
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
fit_reference <- function(reference = NULL, spec, label = "reference") {
  validate_protocol(spec)
  .scalar_string(label, "label")
  if (inherits(reference, "al_imputation")) return(.fit_imputed_reference(reference, spec, label))
  if (spec$method != "cutoff_count" || !is.null(reference)) validate_data(reference, spec)
  if (spec$method == "cutoff_count") return(.fit_cutoffs(reference, spec, label))
  key <- .stratum_key(reference, spec$strata)
  levels <- sort(unique(key))
  if (length(spec$age)) return(.fit_age_reference(reference, spec, label, key, levels))
  rows <- list(); weights <- .reference_weights(reference, spec); weighted <- length(spec$weights) > 0L
  for (level in levels) for (b in spec$biomarkers) {
    x <- .adjusted_values(reference, spec, b)[key == level]; w <- weights[key == level]
    w <- w[!is.na(x)]; x <- x[!is.na(x)]
    if (length(x) < 2L) stop("At least two observed reference values required for ", b$name, .in_stratum(level, spec), call. = FALSE)
    if (spec$method == "quantile_count") {
      p <- .tail_probs(spec, b$direction)
      q <- function(prob) if (weighted) .weighted_quantile(x, w, prob) else unname(stats::quantile(x, prob, type = spec$quantile_type))
      rows[[length(rows) + 1L]] <- .cutoff_row(b, level, spec,
        lower = if ("lower" %in% names(p)) q(p[["lower"]]) else NA_real_,
        upper = if ("upper" %in% names(p)) q(p[["upper"]]) else NA_real_, n = length(x))
    } else {
      tx <- .transform(x, b$transform)
      center <- stats::weighted.mean(tx, w); s <- sqrt(sum(w * (tx - center)^2) / sum(w) * length(tx) / (length(tx) - 1))
      if (!is.finite(s) || s <= 0) stop("Reference SD must be positive for ", b$name, .in_stratum(level, spec), ".", call. = FALSE)
      rows[[length(rows) + 1L]] <- .cutoff_row(b, level, spec, center = center, scale = s, n = length(x))
    }
  }
  .new_reference(spec, do.call(rbind, rows), label, nrow(reference), levels, .strata_values(reference, spec, levels))
}

.fit_cutoffs <- function(reference, spec, label) {
  stratified <- vapply(spec$biomarkers, .stratified_cutoff, logical(1))
  levels <- if (any(stratified)) sort(.cutoff_levels(spec$biomarkers[[which(stratified)[1]]])) else "all"
  key <- if (is.null(reference)) NULL else .stratum_key(reference, spec$strata)
  if (any(stratified) && !is.null(reference) && length(unknown <- setdiff(unique(key), levels)))
    stop("Reference strata without supplied cutoffs: ", paste(unknown, collapse = ", "), call. = FALSE)
  rows <- list()
  for (level in levels) for (b in spec$biomarkers) {
    value <- .cutoff_value(b, level)
    n <- if (is.null(reference)) NA_integer_ else sum(!is.na(reference[[b$name]][key == level]))
    rows[[length(rows) + 1L]] <- .cutoff_row(b, level, spec,
      lower = switch(b$direction, low = value, high = NA_real_, both = value[1]),
      upper = switch(b$direction, low = NA_real_, high = value, both = value[2]), n = n)
  }
  strata_values <- if (any(stratified)) stats::setNames(data.frame(levels), spec$strata) else NULL
  .new_reference(spec, do.call(rbind, rows), label,
    if (is.null(reference)) NA_integer_ else nrow(reference), levels, strata_values)
}

.fit_age_reference <- function(reference, spec, label, key, levels) {
  if (spec$method == "quantile_count") .require_suggested("quantreg", "fit age-conditional quantile references")
  age <- reference[[spec$age]]; weights <- .reference_weights(reference, spec)
  min_n <- 10L * (spec$age_df + 1L)
  rows <- list(); models <- list()
  for (level in levels) for (b in spec$biomarkers) {
    values <- .adjusted_values(reference, spec, b)
    ok <- key == level & !is.na(values)
    where <- paste0(b$name, .in_stratum(level, spec))
    if (sum(ok) < min_n) stop("At least ", min_n, " observed reference values are required for the age model of ", where, ".", call. = FALSE)
    a <- age[ok]; y <- .transform(values[ok], b$transform); w <- weights[ok]
    basis <- splines::ns(a, df = spec$age_df)
    model <- list(biomarker = b$name, stratum = level, transform = b$transform,
      knots = unname(attr(basis, "knots")), boundary = unname(attr(basis, "Boundary.knots")),
      age_range = range(a), coef = list())
    X <- .spline_design(a, model)
    grid <- .spline_design(seq(min(a), max(a), length.out = 101L), model)
    if (spec$method == "quantile_count") {
      p <- .tail_probs(spec, b$direction)
      for (tail in names(p)) model$coef[[tail]] <- .rq_coefficients(X, y, p[[tail]], w)
      if (length(p) == 2L && any(grid %*% model$coef$lower >= grid %*% model$coef$upper))
        stop("Age-specific lower and upper thresholds cross for ", where, ".", call. = FALSE)
    } else {
      fit <- stats::lm.wfit(X, y, w)
      model$coef$center <- unname(fit$coefficients)
      model$coef$scale <- unname(stats::lm.wfit(X, abs(fit$residuals), w)$coefficients) * sqrt(pi / 2)
      if (all(is.finite(model$coef$scale)) && any(grid %*% model$coef$scale <= 0))
        stop("The age model gives non-positive SDs for ", where, ".", call. = FALSE)
    }
    if (any(!is.finite(unlist(model$coef)))) stop("The age model could not be estimated for ", where, ".", call. = FALSE)
    models[[length(models) + 1L]] <- model
    rows[[length(rows) + 1L]] <- .cutoff_row(b, level, spec, n = sum(ok), age_range = model$age_range)
  }
  .new_reference(spec, do.call(rbind, rows), label, nrow(reference), levels,
    .strata_values(reference, spec, levels), models)
}

# Quantile regression coefficients; non-unique solutions are equally valid fits.
.rq_coefficients <- function(X, y, tau, w = NULL) {
  fit <- withCallingHandlers(if (is.null(w) || all(w == w[1])) quantreg::rq.fit(X, y, tau = tau, method = "br")
      else quantreg::rq.wfit(X, y, tau = tau, weights = w, method = "br"),
    warning = function(w) if (grepl("nonunique", conditionMessage(w))) invokeRestart("muffleWarning"))
  unname(fit$coefficients)
}

# Per-row reference parameters for scoring `data`.
.reference_parameters <- function(reference, data) {
  spec <- reference$spec; n <- nrow(data); cut <- reference$cutoffs
  key <- .stratum_key(data, spec$strata)
  if (length(spec$strata) && length(unknown <- setdiff(unique(key), reference$levels)))
    stop("Strata absent from the reference: ", paste(utils::head(unknown, 5L), collapse = ", "), call. = FALSE)
  age <- if (length(spec$age)) data[[spec$age]] else NULL
  lapply(spec$biomarkers, function(b) {
    rows <- which(cut$biomarker == b$name)
    idx <- if (length(spec$strata)) rows[match(key, cut$stratum[rows])] else rep(rows[1L], n)
    column <- function(name) if (name %in% names(cut)) cut[[name]][idx] else rep(NA_real_, n)
    out <- list(lower = column("lower"), upper = column("upper"), center = column("center"),
      scale = column("scale"), outside = rep(FALSE, n), stratum = key)
    if (!is.null(age)) for (model in reference$models) {
      if (model$biomarker != b$name || !any(selected <- key == model$stratum)) next
      a <- age[selected]; G <- .spline_design(a, model)
      for (parameter in names(model$coef)) {
        value <- drop(G %*% model$coef[[parameter]])
        if (parameter %in% c("lower", "upper")) value <- .untransform(value, model$transform)
        out[[parameter]][selected] <- value
      }
      outside <- if (spec$age_extrapolate) rep(FALSE, length(a)) else a < model$age_range[1] | a > model$age_range[2]
      if (!is.null(model$coef$scale)) outside <- outside | out$scale[selected] <= 0
      out$outside[selected] <- outside
    }
    for (parameter in c("lower", "upper", "center", "scale")) out[[parameter]][out$outside] <- NA_real_
    out
  })
}

#' Inspect fitted thresholds
#' @param reference Object from fit_reference().
#' @param at For age-conditional references, a data frame with the age column
#'   (and strata columns) at which to evaluate thresholds. By default, five ages
#'   spanning the reference range are used for every stratum.
#' @return Data frame of thresholds (or reference means and SDs on the
#'   transformed scale for z_mean), units and sample sizes.
#' @family reference calibration
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
reference_cutoffs <- function(reference, at = NULL) {
  if (!inherits(reference, "al_reference")) stop("reference must be an al_reference object.", call. = FALSE)
  spec <- reference$spec
  if (!length(spec$age)) return(reference$cutoffs)
  if (is.null(at)) {
    at <- do.call(rbind, lapply(reference$levels, function(level) {
      d <- reference$cutoffs[if (length(spec$strata)) reference$cutoffs$stratum == level else TRUE, ]
      ages <- seq(max(d$age_min), min(d$age_max), length.out = 5L)
      g <- data.frame(age = ages)
      if (length(spec$strata)) g <- cbind(reference$strata_values[rep(match(level, reference$levels), 5L), , drop = FALSE], g)
      g
    }))
    names(at)[names(at) == "age"] <- spec$age
    rownames(at) <- NULL
  }
  if (!is.data.frame(at) || !all(c(spec$age, spec$strata) %in% names(at)))
    stop("at must be a data frame with the age and strata columns.", call. = FALSE)
  params <- .reference_parameters(reference, at)
  do.call(rbind, lapply(seq_along(spec$biomarkers), function(j) {
    b <- spec$biomarkers[[j]]; p <- params[[j]]
    d <- data.frame(biomarker = b$name, label = b$label, system = b$system, direction = b$direction,
      unit = b$unit)[rep(1L, nrow(at)), ]
    if (length(spec$strata)) d$stratum <- p$stratum
    d$age <- at[[spec$age]]
    d$lower <- p$lower; d$upper <- p$upper
    if (spec$method == "z_mean") { d$center <- p$center; d$scale <- p$scale }
    d$outside_reference <- p$outside
    rownames(d) <- NULL
    d
  }))
}
