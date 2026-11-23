#' Validate a protocol specification
#' @param spec An `al_spec` object.
#' @return TRUE invisibly; malformed specifications cause an error.
#' @family protocol management
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @keywords internal
validate_protocol <- function(spec) {
  if (!inherits(spec, "al_spec")) stop("Expected an define_protocol object.", call. = FALSE)
  for (b in spec$biomarkers) {
    if (!inherits(b, "al_biomarker")) stop("Invalid biomarker object.", call. = FALSE)
    do.call(define_biomarker, b)
  }
  rebuilt <- do.call(define_protocol, unclass(spec))
  if (!identical(rebuilt, spec)) stop("Specification is not canonical; rebuild with define_protocol().", call. = FALSE)
  invisible(TRUE)
}

.assert_reference <- function(x) {
  if (!inherits(x, "al_reference")) stop("Expected an al_reference object.", call. = FALSE)
  validate_protocol(x$spec); spec <- x$spec
  expected <- unname(.biomarker_field(spec, "name"))
  levels <- if (is.null(x$levels)) "all" else x$levels
  if (!is.character(levels) || !length(levels) || anyNA(levels) || anyDuplicated(levels) ||
      (!length(spec$strata) && !identical(levels, "all"))) stop("Invalid reference strata.", call. = FALSE)
  columns <- c("biomarker", "label", "system", "direction", "unit", "lower", "upper", "n_reference",
    if (length(spec$strata)) "stratum", if (spec$method == "z_mean") c("center", "scale"))
  cut <- x$cutoffs
  if (!is.data.frame(cut) || !all(columns %in% names(cut)) ||
      !identical(cut$biomarker, rep(expected, times = length(levels))) ||
      (length(spec$strata) && !identical(cut$stratum, rep(levels, each = length(expected)))))
    stop("Invalid reference threshold table.", call. = FALSE)
  numeric_columns <- intersect(c("lower", "upper", "center", "scale"), names(cut))
  if (!all(vapply(cut[numeric_columns], is.numeric, logical(1)))) stop("Threshold columns must be numeric.", call. = FALSE)
  aged <- length(spec$age) > 0L
  for (i in seq_len(nrow(cut))) {
    d <- cut[i, ]; b <- spec$biomarkers[[match(d$biomarker, expected)]]
    if (!identical(d$unit, b$unit) || !identical(d$direction, b$direction) || !identical(d$system, b$system))
      stop("Reference metadata does not match specification.", call. = FALSE)
    if ((b$direction == "high" && !is.na(d$lower)) || (b$direction == "low" && !is.na(d$upper)))
      stop("Unused threshold contradicts risk direction.", call. = FALSE)
    if (aged) {
      if (!all(is.na(unlist(d[numeric_columns])))) stop("Age-conditional references keep thresholds in models.", call. = FALSE)
      next
    }
    if (spec$method == "z_mean") {
      if (!is.finite(d$center) || !is.finite(d$scale) || d$scale <= 0 || !is.na(d$lower) || !is.na(d$upper))
        stop("Invalid reference mean or SD.", call. = FALSE)
      next
    }
    vals <- c(if (b$direction %in% c("low", "both")) d$lower, if (b$direction %in% c("high", "both")) d$upper)
    if (!is.numeric(vals) || any(!is.finite(vals))) stop("Invalid reference thresholds.", call. = FALSE)
    if (spec$method == "cutoff_count" && !isTRUE(all.equal(vals, .cutoff_value(b, if (is.null(d$stratum)) "all" else d$stratum))))
      stop("Supplied thresholds do not match specification.", call. = FALSE)
    if (b$direction == "both" && d$lower > d$upper) stop("Reversed reference thresholds.", call. = FALSE)
  }
  if (aged) .assert_age_models(x, spec, expected, levels)
  else if (!is.null(x$models)) stop("Only age-conditional references contain models.", call. = FALSE)
  .scalar_string(x$label, "reference label")
  invisible(TRUE)
}

.assert_age_models <- function(x, spec, expected, levels) {
  models <- x$models
  if (!is.list(models) || length(models) != length(expected) * length(levels)) stop("Invalid age models.", call. = FALSE)
  keys <- paste(rep(levels, each = length(expected)), rep(expected, times = length(levels)))
  for (k in seq_along(models)) {
    m <- models[[k]]; b <- spec$biomarkers[[match(m$biomarker, expected)]]
    parameters <- if (spec$method == "z_mean") c("center", "scale") else names(.tail_probs(spec, b$direction))
    if (!is.list(m) || !identical(paste(m$stratum, m$biomarker), keys[k]) || !identical(m$transform, b$transform) ||
        !identical(names(m$coef), parameters) || !is.numeric(m$knots) || any(!is.finite(m$knots)) ||
        length(m$knots) != spec$age_df - 1L || !is.numeric(m$boundary) || length(m$boundary) != 2L ||
        any(!is.finite(m$boundary)) || !is.numeric(m$age_range) || length(m$age_range) != 2L ||
        any(!is.finite(m$age_range)) || m$age_range[1] > m$age_range[2] ||
        !all(vapply(m$coef, function(cf) is.numeric(cf) && length(cf) == spec$age_df + 1L && all(is.finite(cf)), logical(1))))
      stop("Invalid age model for ", m$biomarker, ".", call. = FALSE)
  }
}

# Keep a subset of biomarkers in a fitted reference (panel influence).
.subset_reference <- function(reference, keep) {
  spec <- reference$spec
  names_all <- .biomarker_field(spec, "name")
  kept <- names_all[keep]
  spec$biomarkers <- spec$biomarkers[keep]
  spec$medications <- spec$medications[names(spec$medications) %in% kept]
  reference$spec <- spec
  reference$cutoffs <- reference$cutoffs[reference$cutoffs$biomarker %in% kept, , drop = FALSE]
  rownames(reference$cutoffs) <- NULL
  if (!is.null(reference$models)) reference$models <- Filter(function(m) m$biomarker %in% kept, reference$models)
  # Supplied cutoffs shared by every stratum no longer need stratification.
  if (spec$method == "cutoff_count" && length(spec$strata) &&
      !any(vapply(spec$biomarkers, .stratified_cutoff, logical(1)))) {
    reference$spec$strata <- character()
    totals <- tapply(reference$cutoffs$n_reference, reference$cutoffs$biomarker,
      function(n) if (all(is.na(n))) NA_integer_ else as.integer(sum(n, na.rm = TRUE)))
    first <- reference$cutoffs$stratum == reference$levels[1]
    reference$cutoffs <- reference$cutoffs[first, setdiff(names(reference$cutoffs), "stratum"), drop = FALSE]
    reference$cutoffs$n_reference <- unname(totals[reference$cutoffs$biomarker])
    rownames(reference$cutoffs) <- NULL
    reference$levels <- "all"; reference$strata_values <- NULL
  }
  reference
}

#' Save a versioned protocol passport
#' @param protocol A `define_protocol` or `al_reference`. No participant data are stored.
#' @param file Destination RDS file.
#' @param overwrite Allow replacing an existing file.
#' @return File path invisibly.
#' @family protocol management
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
save_protocol <- function(protocol, file, overwrite = FALSE) {
  if (inherits(protocol, "al_reference")) .assert_reference(protocol) else validate_protocol(protocol)
  .scalar_string(file, "file")
  .check_flag(overwrite, "overwrite")
  if (file.exists(file) && !overwrite) stop("Protocol file already exists.", call. = FALSE)
  saveRDS(list(schema_version = 2L, package_version = .package_version(), protocol = protocol), file)
  invisible(file)
}

#' Load a protocol passport
#' @param file Trusted RDS file saved by save_protocol(). Passports with schema
#'   version 1 (written before strata, age and z-score settings existed) are
#'   upgraded with the new protocol defaults.
#' @return The validated define_protocol or al_reference object.
#' @family protocol management
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
read_protocol <- function(file) {
  .scalar_string(file, "file"); x <- readRDS(file)
  if (!is.list(x) || is.null(x$protocol) || !(identical(x$schema_version, 1L) || identical(x$schema_version, 2L)))
    stop("Unsupported or invalid protocol schema.", call. = FALSE)
  protocol <- x$protocol
  missing_fields <- function(spec) inherits(spec, "al_spec") && !all(names(formals(define_protocol)) %in% names(spec))
  if (identical(x$schema_version, 2L) && inherits(protocol, "al_reference") && missing_fields(protocol$spec)) protocol$spec <- .upgrade_spec(protocol$spec)
  if (identical(x$schema_version, 2L) && missing_fields(protocol)) protocol <- .upgrade_spec(protocol)
  if (identical(x$schema_version, 1L)) {
    if (inherits(protocol, "al_reference")) {
      protocol$spec <- .upgrade_spec(protocol$spec)
      protocol <- .new_reference(protocol$spec, protocol$cutoffs, protocol$label, protocol$n_reference, "all")
    } else if (inherits(protocol, "al_spec")) protocol <- .upgrade_spec(protocol)
  }
  if (inherits(protocol, "al_reference")) .assert_reference(protocol) else validate_protocol(protocol)
  protocol
}

.upgrade_spec <- function(spec) {
  if (!inherits(spec, "al_spec")) stop("Unsupported or invalid protocol schema.", call. = FALSE)
  old <- unclass(spec)
  old$biomarkers <- lapply(old$biomarkers, function(b) do.call(define_biomarker, unclass(b)))
  do.call(define_protocol, old)
}

.named_list <- function(x, min_length = 1L) {
  if (!is.list(x) || length(x) < min_length || is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))
    stop("Supply a list with unique, non-empty names.", call. = FALSE)
}

#' Check comparability of fitted references
#' @param references Named list of at least two al_reference objects.
#' @return Pairwise data frame of algorithm and threshold compatibility. This
#'   does not establish population exchangeability or clinical validity.
#' @family protocol management
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
compare_references <- function(references) {
  .named_list(references, 2L); invisible(lapply(references, .assert_reference))
  pairs <- utils::combn(seq_along(references), 2)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    a <- references[[pairs[1, i]]]; b <- references[[pairs[2, i]]]
    config <- function(r) {
      s <- r$spec
      meds <- s$medications
      if (length(meds)) meds <- meds[order(names(meds))] else meds <- character()
      list(method = s$method, missing = s$missing, min_coverage = s$min_coverage,
        inclusive = s$inclusive, medications = meds, balance = s$balance, bilateral = s$bilateral,
        strata = s$strata, age = s$age, age_df = s$age_df, age_extrapolate = s$age_extrapolate, z_cap = s$z_cap,
        medication_rule = s$medication_rule, medication_effects = s$medication_effects, weights = s$weights)
    }
    panel <- function(r) {
      d <- do.call(rbind, lapply(r$spec$biomarkers, function(m) data.frame(biomarker = m$name,
        system = m$system, direction = m$direction, unit = m$unit, transform = m$transform)))
      d <- d[order(d$biomarker), ]; rownames(d) <- NULL; d
    }
    cut <- function(r) {
      columns <- intersect(c("biomarker", "stratum", "lower", "upper", "center", "scale"), names(r$cutoffs))
      d <- r$cutoffs[do.call(order, unname(as.list(r$cutoffs[intersect(c("biomarker", "stratum"), columns)]))), columns]
      rownames(d) <- NULL
      models <- r$models
      if (length(models)) models <- models[order(vapply(models, function(m) paste(m$biomarker, m$stratum), character(1)))]
      list(d, models)
    }
    panel_ok <- identical(panel(a), panel(b))
    algorithm_ok <- panel_ok && identical(config(a), config(b))
    thresholds_ok <- panel_ok && identical(cut(a), cut(b))
    reasons <- c(if (!panel_ok) "panel/units/directions/systems/transforms differ",
                 if (!algorithm_ok) "scoring rules differ",
                 if (!thresholds_ok) "thresholds differ")
    data.frame(reference_a = names(references)[pairs[1, i]], reference_b = names(references)[pairs[2, i]],
      panel_equal = panel_ok, algorithm_equal = algorithm_ok, thresholds_equal = thresholds_ok,
      compatible = algorithm_ok && thresholds_ok,
      reason = if (length(reasons)) paste(reasons, collapse = "; ") else "same scoring rules and thresholds")
  }))
}
