#' Score data against a fixed reference
#'
#' Applies fixed reference parameters without refitting. Each biomarker
#' contributes a risk flag (count methods) or a risk-oriented z-score
#' (`z_mean`; high-risk direction positive). Contributions are aggregated per
#' the protocol: counts, proportions among observed biomarkers, mean z-scores,
#' or system-balanced means in which each system weighs equally.
#'
#' @param data Biomarker data frame, in the units declared by the protocol.
#' @param reference Fitted reference from fit_reference().
#' @param id Optional participant identifier column to retain in the scores table.
#' @param group Optional column to retain for descriptive group graphics.
#' @return An `al_result` with scores, components, systems, audit metadata and plots.
#'   The score scale is recorded in `audit$scale`: `count`, `proportion`,
#'   `system_proportion`, `z` or `system_z`. Status `age_outside_reference`
#'   marks participants outside the age range of an age-conditional reference.
#' @family load scoring
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
score_load <- function(data, reference, id = NULL, group = NULL) {
  .assert_reference(reference)
  spec <- reference$spec; validate_data(data, spec)
  for (col in list(id, group)) if (!is.null(col)) {
    .scalar_string(col, "metadata column")
    if (!col %in% names(data)) stop("Missing metadata column: ", col, call. = FALSE)
    if (!is.atomic(data[[col]]) || !is.null(dim(data[[col]]))) stop("Metadata must be a vector.", call. = FALSE)
  }
  for (column in .invalid_medication_columns(data, spec$medications))
    stop("Medication columns must be complete logical vectors: ", column, call. = FALSE)
  ids <- if (is.null(id)) seq_len(nrow(data)) else data[[id]]
  if (anyNA(ids) || anyDuplicated(ids)) stop("Identifiers must be unique and non-missing.", call. = FALSE)
  params <- .reference_parameters(reference, data)
  count <- .is_count_method(spec)
  components <- lapply(seq_along(spec$biomarkers), function(j) {
    b <- spec$biomarkers[[j]]; p <- params[[j]]; raw <- data[[b$name]]; x <- .adjusted_values(data, spec, b); n <- length(x)
    medication <- adjusted <- rep(FALSE, n); risk <- rep(NA_integer_, n)
    if (count) {
      flag <- .risk_flag(x, p, b$direction, spec$inclusive)
      if (b$name %in% names(spec$medications)) medication <- data[[spec$medications[[b$name]]]]
      if (spec$medication_rule == "flag") {
        adjusted <- medication & !is.na(x) & !flag
        adjusted[is.na(adjusted)] <- FALSE
        flag <- flag | medication
      } else adjusted <- medication & !is.na(x)
      flag[is.na(x) | p$outside] <- NA
      adjusted[is.na(flag)] <- FALSE
      risk <- as.integer(flag); contribution <- as.numeric(risk)
    } else {
      contribution <- (.transform(x, b$transform) - p$center) / p$scale
      if (b$direction == "low") contribution <- -contribution
      contribution <- pmin(pmax(contribution, -spec$z_cap), spec$z_cap)
      contribution[p$outside] <- NA_real_
      if (b$name %in% names(spec$medications)) {
        medication <- data[[spec$medications[[b$name]]]]; adjusted <- medication & !is.na(x)
      }
    }
    d <- data.frame(row = seq_len(n), id = ids, biomarker = b$name,
      label = b$label, system = b$system, unit = b$unit, direction = b$direction,
      value = raw, lower = p$lower, upper = p$upper, risk = risk,
      medication_flag = medication, medication_adjusted = adjusted)
    if (spec$medication_rule == "add") d$value_adjusted <- x
    if (!count) { d$center <- p$center; d$scale <- p$scale; d$z <- contribution }
    d$contribution <- contribution
    if (length(spec$strata)) d$stratum <- p$stratum
    d
  })
  mat <- do.call(cbind, lapply(components, `[[`, "contribution"))
  observed <- rowSums(!is.na(mat)); total <- ncol(mat); coverage <- observed / total
  enough <- if (spec$missing == "complete") observed == total else coverage >= spec$min_coverage & observed > 0
  system_of <- .biomarker_field(spec, "system")
  system_levels <- unique(system_of)
  system_mat <- matrix(vapply(system_levels, function(s) {
    sub <- mat[, system_of == s, drop = FALSE]; n <- rowSums(!is.na(sub))
    value <- rowSums(sub, na.rm = TRUE) / pmax(n, 1); value[n == 0] <- NA_real_; value
  }, numeric(nrow(mat))), nrow = nrow(mat))
  score <- if (spec$balance == "system") rowMeans(system_mat)
    else if (count && spec$missing == "complete") rowSums(mat)
    else rowSums(mat, na.rm = TRUE) / pmax(observed, 1)
  valid <- enough & !is.na(score)
  score[!valid] <- NA_real_
  outside <- Reduce(`|`, lapply(params, `[[`, "outside"))
  scores <- data.frame(row = seq_len(nrow(data)), id = ids, score = score,
    n_observed = observed, n_total = total, coverage = coverage,
    status = ifelse(valid, "ok", ifelse(outside, "age_outside_reference", "insufficient_coverage")))
  if (!is.null(group)) scores$group <- data[[group]]
  if (length(spec$strata)) scores$stratum <- params[[1]]$stratum
  if (length(spec$age)) scores$age <- data[[spec$age]]
  cmp <- do.call(rbind, components); rownames(cmp) <- NULL
  systems <- do.call(rbind, lapply(seq_along(system_levels), function(k) {
    s <- system_levels[k]; n <- rowSums(!is.na(mat[, system_of == s, drop = FALSE]))
    data.frame(row = scores$row, id = ids, system = s, score = system_mat[, k],
      risk_proportion = if (count) system_mat[, k] else NA_real_,
      n_observed = n, n_total = sum(system_of == s))
  }))
  structure(list(scores = scores, components = cmp, systems = systems,
    reference = reference, audit = list(method = spec$method, missing = spec$missing,
      balance = spec$balance, inclusive = spec$inclusive, min_coverage = spec$min_coverage,
      scale = .score_scale(spec), strata = spec$strata, age = spec$age,
      package_version = .package_version(), declared_units = attr(data, "al_units"),
      conversions = attr(data, "al_conversions"), units = "Unannotated columns assumed to use protocol units"),
    plots = NULL), class = "al_result")
}

#' Calculate allostatic load with optional graphics
#' @param data Biomarker data frame.
#' @param spec Protocol; required when reference is not supplied.
#' @param reference Optional fitted external reference.
#' @param graphics No plots by default (none); auto creates the optional plot set.
#' @param id,group Optional identifier and grouping columns.
#' @param locale Plot label language: English (en).
#' @return An `al_result`. When graphics is auto, plots contains editable ggplots.
#' @examples
#' source(system.file("examples", "panel.R", package = "allostaticR"))
#' result <- allostatic_load(data, spec, id = "id", group = "group")
#' summary(result)
#' plot_data(result)
#' plot_load(result, type = "boxplot")
#' @family load scoring
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
allostatic_load <- function(data, spec = NULL, reference = NULL,
                            graphics = c("none", "auto"), id = NULL,
                            group = NULL, locale = "en") {
  graphics <- match.arg(graphics); locale <- match.arg(locale, "en")
  if (!is.null(reference) && !is.null(spec)) stop("Supply spec or reference, not both.", call. = FALSE)
  if (is.null(reference)) {
    if (is.null(spec)) stop("Supply spec or reference.", call. = FALSE)
    reference <- fit_reference(data, spec, label = if (spec$method == "cutoff_count") "supplied_cutoffs" else "study_sample")
  }
  result <- score_load(data, reference, id, group)
  if (graphics == "auto") result$plots <- .plot_set(result, locale = locale)
  result
}
