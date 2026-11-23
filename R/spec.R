.scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop(name, " must be a non-empty string.", call. = FALSE)
}

# Stratified cutoffs: named numbers (high/low) or a named list of ordered pairs (both).
.stratified_cutoff <- function(b) {
  !is.null(b$cutoff) && if (b$direction == "both") is.list(b$cutoff) else !is.null(names(b$cutoff))
}
.cutoff_levels <- function(b) if (.stratified_cutoff(b)) names(b$cutoff) else character()
.cutoff_value <- function(b, level) {
  if (!.stratified_cutoff(b)) return(unname(unlist(b$cutoff)))
  unname(if (is.list(b$cutoff)) b$cutoff[[level]] else b$cutoff[[level]])
}

.check_cutoff <- function(cutoff, direction, name) {
  pair <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x)) && x[1] < x[2]
  ok <- if (direction == "both") {
    if (is.list(cutoff)) length(cutoff) > 0L && .valid_names(names(cutoff)) && all(vapply(cutoff, pair, logical(1)))
    else pair(cutoff)
  } else is.numeric(cutoff) && length(cutoff) >= 1L && all(is.finite(cutoff)) &&
    (if (is.null(names(cutoff))) length(cutoff) == 1L else .valid_names(names(cutoff)))
  if (!ok) stop("Invalid cutoff for biomarker ", name, ".", call. = FALSE)
}

#' Define a biomarker
#' @param name Column name.
#' @param system Physiological system.
#' @param direction Risk direction: high, low, or both.
#' @param unit Declared measurement unit. No conversion is performed.
#' @param cutoff Numeric cutoff, or two ordered limits for bilateral risk. For
#'   stratum-specific supplied cutoffs (for example sex-specific HDL limits),
#'   use a named numeric vector (`c(Female = 50, Male = 40)`) or, for bilateral
#'   risk, a named list of ordered pairs. Names are levels of the single
#'   `strata` column declared in [define_protocol()].
#' @param label Display label.
#' @param transform Scale used for z-scores and age models: `"identity"` or
#'   `"log"` (positive values only). Thresholds are always reported in the
#'   declared unit; quantile flags are unaffected by monotone transforms.
#' @return An `al_biomarker` specification.
#' @family protocol definition
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
define_biomarker <- function(name, system, direction = c("high", "low", "both"),
                         unit, cutoff = NULL, label = name, transform = c("identity", "log")) {
  .scalar_string(name, "name"); .scalar_string(system, "system")
  .scalar_string(unit, "unit"); .scalar_string(label, "label")
  direction <- match.arg(direction); transform <- match.arg(transform)
  if (!is.null(cutoff)) .check_cutoff(cutoff, direction, name)
  structure(list(name = name, system = system, direction = direction,
                 unit = unit, cutoff = cutoff, label = label, transform = transform), class = "al_biomarker")
}

#' Define an allostatic load protocol
#'
#' A protocol fixes every scoring decision so that it can be fitted once and
#' reused: how biomarker contributions are formed, how they are aggregated,
#' which reference covariates condition the thresholds, and how missing
#' measurements are handled.
#'
#' @param biomarkers List of objects from define_biomarker().
#' @param method Contribution of each biomarker: `"quantile_count"` (risk flag
#'   from reference quantiles), `"cutoff_count"` (risk flag from supplied
#'   cutoffs) or `"z_mean"` (risk-oriented z-score against the reference mean
#'   and SD; high and low directions only).
#' @param risk_quantile Upper-tail probability; the lower tail uses its complement.
#' @param missing Complete panels or proportion among observed biomarkers.
#' @param min_coverage Minimum observed fraction for proportion scores.
#' @param inclusive Whether equality with a cutoff denotes risk.
#' @param quantile_type Quantile algorithm passed to stats::quantile().
#' @param protocol_id,protocol_version User-defined protocol identifier and version.
#' @param sources Character vector of protocol citations or URLs.
#' @param medications Named character vector mapping biomarkers to complete logical
#'   medication columns. Medication flags may not contain NA; missing
#'   measurements remain NA.
#' @param medication_rule How treatment enters the score. `"flag"`: a TRUE
#'   flag marks a measured biomarker as at risk (count methods only). `"add"`:
#'   a fixed treatment effect is added to the measured value of treated
#'   participants before thresholds or z-scores are applied, in the reference
#'   and in scored data (for example +10 mmHg systolic blood pressure for
#'   antihypertensive treatment; Tobin et al. 2005).
#' @param medication_effects For `medication_rule = "add"`, a named numeric
#'   vector of effects in declared units, one per medication biomarker. Effects
#'   must come from the literature; none are supplied by the package.
#' @param weights Optional column with positive sampling weights used when
#'   fitting the reference (weighted quantiles, means, SDs and age models), as
#'   in NHANES-type surveys. Scored data do not need it.
#' @param balance `"biomarker"` gives every biomarker equal weight;
#'   `"system"` averages within each physiological system first, so every
#'   system has equal weight regardless of how many biomarkers it contains.
#' @param bilateral For quantile references of bilateral biomarkers: `"full"`
#'   uses both tails at the risk quantile (about 50% flagged at 0.75);
#'   `"split"` divides the tail probability between both tails (about 25%).
#' @param strata Columns defining reference strata (for example `"sex"`).
#'   Thresholds or means are fitted within each stratum and scored data must
#'   contain the same complete columns; unknown strata are rejected.
#' @param age Optional numeric column for age-conditional references: quantile
#'   regression on a natural spline of age for `quantile_count` (requires the
#'   quantreg package) and the absolute-residual mean/SD model of Altman (1993)
#'   for `z_mean`.
#' @param age_df Degrees of freedom of the natural spline of age.
#' @param age_extrapolate Whether to evaluate age models outside the observed
#'   reference age range. When FALSE, such participants receive NA
#'   contributions and the status `age_outside_reference`.
#' @param z_cap Optional absolute cap for z-scores (`Inf` for none).
#' @return An `al_spec` protocol object.
#' @family protocol definition
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
define_protocol <- function(biomarkers, method = c("quantile_count", "cutoff_count", "z_mean"),
                    risk_quantile = 0.75, missing = c("complete", "proportion"),
                    min_coverage = 1, inclusive = TRUE, quantile_type = 7L,
                    protocol_id = "custom", protocol_version = "1",
                    sources = character(), medications = character(),
                    balance = c("biomarker", "system"), bilateral = c("full", "split"),
                    strata = character(), age = character(), age_df = 3L,
                    age_extrapolate = FALSE, z_cap = Inf, medication_rule = c("flag", "add"),
                    medication_effects = numeric(), weights = character()) {
  method <- match.arg(method); missing <- match.arg(missing)
  balance <- match.arg(balance); bilateral <- match.arg(bilateral); medication_rule <- match.arg(medication_rule)
  if (!is.list(biomarkers) || !length(biomarkers) ||
      !all(vapply(biomarkers, inherits, logical(1), "al_biomarker")))
    stop("biomarkers must be a non-empty list of define_biomarker objects.", call. = FALSE)
  nm <- vapply(biomarkers, `[[`, character(1), "name")
  if (anyDuplicated(nm)) stop("Biomarker names must be unique.", call. = FALSE)
  if (!is.numeric(risk_quantile) || length(risk_quantile) != 1L ||
      !is.finite(risk_quantile) || risk_quantile <= 0.5 || risk_quantile >= 1)
    stop("risk_quantile must be between 0.5 and 1.", call. = FALSE)
  if (!is.numeric(min_coverage) || length(min_coverage) != 1L ||
      !is.finite(min_coverage) || min_coverage <= 0 || min_coverage > 1)
    stop("min_coverage must be in (0, 1].", call. = FALSE)
  if (!is.logical(inclusive) || length(inclusive) != 1L || is.na(inclusive))
    stop("inclusive must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(quantile_type) || length(quantile_type) != 1L ||
      is.na(quantile_type) || !quantile_type %in% 1:9)
    stop("quantile_type must be an integer from 1 to 9.", call. = FALSE)
  if (method == "cutoff_count" && any(vapply(biomarkers, function(b) is.null(b$cutoff), logical(1))))
    stop("Every biomarker needs a cutoff for cutoff_count.", call. = FALSE)
  .scalar_string(protocol_id, "protocol_id"); .scalar_string(protocol_version, "protocol_version")
  if (!is.character(sources) || anyNA(sources)) stop("sources must be a character vector without NA.", call. = FALSE)
  if (!is.character(medications) || anyNA(medications) ||
      (length(medications) && (is.null(names(medications)) || anyNA(names(medications)) ||
       anyDuplicated(names(medications)) || any(!names(medications) %in% nm) || any(!nzchar(medications)))))
    stop("medications must map unique biomarker names to column names.", call. = FALSE)
  if (method == "z_mean") {
    if (any(vapply(biomarkers, `[[`, character(1), "direction") == "both"))
      stop("z_mean supports high and low risk directions only.", call. = FALSE)
    if (length(medications) && medication_rule == "flag")
      stop("Medication flags require a count method; use medication_rule = 'add' with z_mean.", call. = FALSE)
  }
  if (medication_rule == "add" && !length(medications) && !length(medication_effects)) {
    medication_effects <- numeric()
  } else if (medication_rule == "add") {
    if (!is.numeric(medication_effects) || length(medication_effects) != length(medications) ||
        !.valid_names(names(medication_effects)) || !setequal(names(medication_effects), names(medications)) ||
        any(!is.finite(medication_effects)))
      stop("medication_effects must give one finite effect per medication biomarker.", call. = FALSE)
    medication_effects <- medication_effects[names(medications)]
  } else if (length(medication_effects)) stop("medication_effects require medication_rule = 'add'.", call. = FALSE)
  if (!is.character(weights) || length(weights) > 1L || anyNA(weights) || any(!nzchar(weights)))
    stop("weights must be a single column name or character().", call. = FALSE)
  if (!is.character(strata) || anyNA(strata) || any(!nzchar(strata)) || anyDuplicated(strata))
    stop("strata must be a character vector of unique column names.", call. = FALSE)
  if (!is.character(age) || length(age) > 1L || anyNA(age) || any(!nzchar(age)))
    stop("age must be a single column name or character().", call. = FALSE)
  if (any(c(strata, age, weights) %in% c(nm, unname(medications))) || any(weights %in% c(strata, age)))
    stop("strata and age columns cannot be biomarker or medication columns.", call. = FALSE)
  if (length(age) && age %in% strata) stop("age cannot also be a stratum.", call. = FALSE)
  if (!is.numeric(age_df) || length(age_df) != 1L || is.na(age_df) || age_df != floor(age_df) ||
      age_df < 1 || age_df > 10) stop("age_df must be an integer from 1 to 10.", call. = FALSE)
  if (!is.logical(age_extrapolate) || length(age_extrapolate) != 1L || is.na(age_extrapolate))
    stop("age_extrapolate must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(z_cap) || length(z_cap) != 1L || is.na(z_cap) || z_cap <= 0)
    stop("z_cap must be a positive number or Inf.", call. = FALSE)
  if (method == "cutoff_count") {
    if (length(age)) stop("Age-conditional references require quantile_count or z_mean.", call. = FALSE)
    stratified <- vapply(biomarkers, .stratified_cutoff, logical(1))
    if (any(stratified)) {
      if (length(strata) != 1L) stop("Stratified cutoffs require exactly one strata column.", call. = FALSE)
      levels <- lapply(biomarkers[stratified], .cutoff_levels)
      if (!all(vapply(levels, setequal, logical(1), levels[[1]])))
        stop("Stratified cutoffs must share the same stratum levels.", call. = FALSE)
    } else if (length(strata)) stop("strata with cutoff_count require stratified cutoffs.", call. = FALSE)
  }
  structure(list(biomarkers = biomarkers, method = method,
    risk_quantile = risk_quantile, missing = missing, min_coverage = min_coverage,
    inclusive = inclusive, quantile_type = quantile_type, protocol_id = protocol_id,
    protocol_version = protocol_version, sources = sources, medications = medications,
    balance = balance, bilateral = bilateral, strata = strata, age = age, age_df = age_df,
    age_extrapolate = age_extrapolate, z_cap = z_cap, medication_rule = medication_rule,
    medication_effects = medication_effects, weights = weights), class = "al_spec")
}

#' Validate biomarker data
#' @param data Data frame containing biomarker columns and, when the protocol
#'   declares them, complete strata and age columns.
#' @param spec Protocol object (`al_spec`) created by [define_protocol()].
#' @return TRUE invisibly, or an informative error.
#' @family protocol definition
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @keywords internal
validate_data <- function(data, spec) {
  validate_protocol(spec)
  .check_data_frame(data)
  units <- attr(data, "al_units"); .check_units_vector(units)
  for (b in spec$biomarkers) {
    if (b$name %in% names(units) && units[[b$name]] != b$unit)
      stop("Declared unit mismatch for ", b$name, call. = FALSE)
    if (!b$name %in% names(data)) stop("Missing biomarker column: ", b$name, call. = FALSE)
    x <- data[[b$name]]
    if (!is.numeric(x) || !is.null(dim(x))) stop(b$name, " must be a numeric vector.", call. = FALSE)
    if (any(is.infinite(x)) || any(is.nan(x))) stop(b$name, " contains non-finite values other than NA.", call. = FALSE)
    if (b$transform == "log" && any(x <= 0, na.rm = TRUE))
      stop(b$name, " must be positive for the log transform.", call. = FALSE)
  }
  for (column in spec$strata) {
    if (!column %in% names(data)) stop("Missing strata column: ", column, call. = FALSE)
    x <- data[[column]]
    if (!is.atomic(x) || !is.null(dim(x)) || anyNA(x))
      stop("Strata column ", column, " must be a complete vector.", call. = FALSE)
  }
  if (length(spec$age)) {
    if (!spec$age %in% names(data)) stop("Missing age column: ", spec$age, call. = FALSE)
    x <- data[[spec$age]]
    if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x)))
      stop("Age column ", spec$age, " must be complete and finite.", call. = FALSE)
  }
  invisible(TRUE)
}
