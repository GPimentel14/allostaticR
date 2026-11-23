.templates <- list(
  seeman1997 = list(
    citation = paste("Seeman TE, Singer BH, Rowe JW, Horwitz RI, McEwen BS (1997). Price of adaptation:",
      "allostatic load and its health consequences. MacArthur studies of successful aging.",
      "Archives of Internal Medicine 157, 2259-2268. doi:10.1001/archinte.1997.00440400111013"),
    note = "Original rule: count of biomarkers in the high-risk quartile (lowest quartile for HDL and DHEA-S).",
    biomarkers = data.frame(
      name = c("sbp", "dbp", "waist_hip_ratio", "total_hdl_ratio", "hdl", "hba1c", "urinary_cortisol",
        "urinary_norepinephrine", "urinary_epinephrine", "dhea_s"),
      label = c("Systolic blood pressure", "Diastolic blood pressure", "Waist-hip ratio", "Total/HDL cholesterol ratio",
        "HDL cholesterol", "Glycated hemoglobin", "Urinary cortisol (12 h)", "Urinary norepinephrine (12 h)",
        "Urinary epinephrine (12 h)", "DHEA-S"),
      system = c("cardiovascular", "cardiovascular", "metabolic", "metabolic", "metabolic", "metabolic",
        "neuroendocrine", "sympathetic", "sympathetic", "neuroendocrine"),
      direction = c("high", "high", "high", "high", "low", "high", "high", "high", "high", "low"),
      unit = c("mmHg", "mmHg", "ratio", "ratio", "mg/dL", "%", "ug/g creatinine", "ug/g creatinine",
        "ug/g creatinine", "ug/dL"),
      transform = c("identity", "identity", "identity", "identity", "identity", "identity", "log", "log", "log", "log"))))

#' Build or list published protocol templates
#'
#' Reproduces the biomarker panel, systems and risk directions of a published
#' allostatic load operationalization, with its citation recorded in
#' `sources`. Calling `protocol_template()` with `name = NULL` returns a data
#' frame listing all available templates. Templates fix the panel, not the
#' thresholds: fit them to your own reference sample (the original rule used
#' sample quartiles). Units are typical declarations only; replace them with
#' the units of your data.
#'
#' @param name Template name, or `NULL` (default) to list available templates.
#' @param columns Optional named character vector mapping template biomarker
#'   names to your column names, for example `c(dhea_s = "dheas_ugdl")`.
#' @param units Optional named character vector overriding declared units.
#' @param drop Optional template biomarkers to omit (for example when a
#'   measurement is unavailable); omissions are recorded in `sources`.
#' @param ... Further arguments for [define_protocol()], such as `method`,
#'   `strata` or `missing`.
#' @return A data frame of available template biomarkers when `name` is `NULL`,
#'   otherwise an `al_spec` protocol.
#' @examples
#' protocol_template()[, c("template", "name", "direction", "unit")]
#' spec <- protocol_template("seeman1997", drop = c("urinary_norepinephrine", "urinary_epinephrine"))
#' length(spec$biomarkers)
#' @family protocol definition
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
protocol_template <- function(name = NULL, columns = character(), units = character(), drop = character(), ...) {
  if (is.null(name)) {
    return(do.call(rbind, lapply(names(.templates), function(t) cbind(template = t, .templates[[t]]$biomarkers,
      citation = .templates[[t]]$citation))))
  }
  if (!is.character(name) || length(name) != 1L || !name %in% names(.templates))
    stop("Unknown template; available: ", paste(names(.templates), collapse = ", "), call. = FALSE)
  t <- .templates[[name]]; d <- t$biomarkers
  for (arg in list(columns, units)) if (length(arg) && (!is.character(arg) || !.valid_names(names(arg)) || !all(names(arg) %in% d$name)))
    stop("columns and units must be named by template biomarkers.", call. = FALSE)
  if (!is.character(drop) || !all(drop %in% d$name)) stop("drop must name template biomarkers.", call. = FALSE)
  d <- d[!d$name %in% drop, , drop = FALSE]
  if (!nrow(d)) stop("No biomarkers left.", call. = FALSE)
  d$column <- ifelse(d$name %in% names(columns), columns[d$name], d$name)
  d$unit <- ifelse(d$name %in% names(units), units[d$name], d$unit)
  biomarkers <- lapply(seq_len(nrow(d)), function(i) define_biomarker(unname(d$column[i]), d$system[i], d$direction[i],
    unname(d$unit[i]), label = d$label[i], transform = d$transform[i]))
  sources <- c(t$citation, t$note, if (length(drop)) paste("Omitted from template:", paste(drop, collapse = ", ")))
  define_protocol(biomarkers, protocol_id = name, sources = sources, ...)
}
