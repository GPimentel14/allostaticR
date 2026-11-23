#' Audit biomarker data quality
#' @param data Data frame to audit; invalid columns are reported rather than repaired.
#' @inheritParams validate_data
#' @param limits Optional named list of two finite increasing plausibility bounds
#'   per biomarker. These are user-defined limits, not universal clinical ranges.
#' @param units Optional named character vector declaring actual column units.
#' @param id Optional column to check for missing and duplicate identifiers.
#' @return An `al_quality` object with per-biomarker quality statistics, explicit
#'   issues and a missingness graphic. No data are changed or imputed.
#' @family data quality
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
check_data_quality <- function(data, spec, limits = list(), units = NULL, id = NULL) {
  validate_protocol(spec)
  .check_data_frame(data)
  markers <- .biomarker_field(spec, "name")
  if (!is.list(limits) || (length(limits) && (is.null(names(limits)) || anyNA(names(limits)) || anyDuplicated(names(limits)) || any(!names(limits) %in% markers))))
    stop("limits must be a named biomarker list.", call. = FALSE)
  for (x in limits) if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || x[1] >= x[2])
    stop("Each plausibility limit needs two finite increasing bounds.", call. = FALSE)
  if (is.null(units)) units <- attr(data, "al_units")
  .check_units_vector(units)
  stats <- .per_biomarker(spec, function(b) {
    present <- b$name %in% names(data); x <- if (present) data[[b$name]] else rep(NA_real_, nrow(data))
    numeric <- is.numeric(x) && is.null(dim(x))
    finite <- if (numeric) is.finite(x) else rep(FALSE, nrow(data))
    bounds <- limits[[b$name]]
    outside <- if (numeric && !is.null(bounds)) sum(finite & (x < bounds[1] | x > bounds[2])) else NA_integer_
    declared <- if (!is.null(units) && b$name %in% names(units)) units[[b$name]] else NA_character_
    data.frame(biomarker = b$name, present = present, numeric = numeric,
      n = nrow(data), n_missing = if (present) sum(is.na(x)) else nrow(data),
      n_nonfinite = if (numeric) sum(is.infinite(x) | is.nan(x)) else NA_integer_,
      n_outside_limits = outside, expected_unit = b$unit, declared_unit = declared,
      unit_status = if (is.na(declared)) "undeclared" else if (identical(declared, b$unit)) "ok" else "mismatch")
  })
  issues <- c(if (any(!stats$present)) "missing_columns", if (any(!stats$numeric)) "non_numeric_columns",
    if (any(stats$n_nonfinite > 0, na.rm = TRUE)) "nonfinite_values",
    if (any(stats$n_outside_limits > 0, na.rm = TRUE)) "outside_user_limits",
    if (any(stats$unit_status == "mismatch")) "unit_mismatch",
    if (any(stats$unit_status == "undeclared")) "units_not_verified")
  for (column in .invalid_medication_columns(data, spec$medications))
    issues <- c(issues, paste0("invalid_medication_column:", column))
  if (!is.null(id)) {
    .scalar_string(id, "id")
    if (!id %in% names(data)) stop("Missing identifier column.", call. = FALSE)
    if (anyNA(data[[id]])) issues <- c(issues, "missing_identifiers")
    if (anyDuplicated(data[[id]])) issues <- c(issues, "duplicate_identifiers")
  }
  stats$missing_fraction <- stats$n_missing / stats$n
  p <- .finish_plot(
    ggplot2::ggplot(stats, ggplot2::aes(x = .data$missing_fraction, y = .data$biomarker)) +
      ggplot2::geom_col(fill = .al_palette$groups[["ink"]], width = 0.65) +
      ggplot2::scale_x_continuous(limits = c(0, 1), labels = .percent),
    title = "Biomarker missingness", x = "Missing observations (%)",
    y = "Biomarker", caption = "Source: supplied data. Absent columns counted as fully unavailable; no imputation.")
  structure(list(summary = stats, issues = issues, plots = list(missingness = p)), class = "al_quality")
}

.check_units_vector <- function(units) {
  if (!is.null(units) && (!is.character(units) || is.null(names(units)) ||
      anyNA(units) || anyNA(names(units)) || any(!nzchar(units)) ||
      any(!nzchar(names(units))) || anyDuplicated(names(units))))
    stop("units must be a named character vector of non-empty units.", call. = FALSE)
}

#' Convert a column with an explicit unit transformation
#' @param data Data frame.
#' @param column Numeric biomarker column.
#' @param from,to Source and target units.
#' @param factor Positive finite multiplier supplied by the user.
#' @param offset Finite additive offset; converted = original * factor + offset.
#' @return Converted data frame with al_units and al_conversions attributes.
#' @details The caller must verify the conversion factor for the analyte.
#'   This function does not guess factors from unit names. Declared source-unit
#'   conflicts are rejected. Missing values remain missing.
#' @family data quality
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
convert_units <- function(data, column, from, to, factor, offset = 0) {
  .scalar_string(column, "column"); .scalar_string(from, "from"); .scalar_string(to, "to")
  if (!is.data.frame(data) || !column %in% names(data) || !is.numeric(data[[column]]) ||
      !is.null(dim(data[[column]])) || anyDuplicated(names(data))) stop("Expected a unique numeric column.", call. = FALSE)
  if (!is.numeric(factor) || length(factor) != 1L || !is.finite(factor) || factor <= 0 ||
      !is.numeric(offset) || length(offset) != 1L || !is.finite(offset)) stop("Invalid conversion factor/offset.", call. = FALSE)
  if (any(is.infinite(data[[column]]) | is.nan(data[[column]]))) stop("Cannot convert non-finite values other than NA.", call. = FALSE)
  units <- attr(data, "al_units"); .check_units_vector(units)
  if (column %in% names(units) && units[[column]] != from) stop("Declared source unit mismatch.", call. = FALSE)
  converted <- data[[column]] * factor + offset
  if (any(is.infinite(converted))) stop("Conversion overflow.", call. = FALSE)
  data[[column]] <- converted; units[column] <- to
  attr(data, "al_units") <- units
  attr(data, "al_conversions") <- c(attr(data, "al_conversions"), list(
    list(column = column, from = from, to = to, factor = factor, offset = offset)))
  data
}
