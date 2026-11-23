.assert_result <- function(x) {
  if (!inherits(x, "al_result")) stop("Expected an al_result object.", call. = FALSE)
}

#' @export
print.al_result <- function(x, ...) {
  cat("Allostatic load:", x$audit$method, "| scale:", x$audit$scale,
    if (length(x$audit$strata)) paste("| strata:", paste(x$audit$strata, collapse = ", ")),
    if (length(x$audit$age)) paste("| age-conditional:", x$audit$age), "\n")
  cat(sum(!is.na(x$scores$score)), "of", nrow(x$scores), "records scored | reference:", x$reference$label, "\n")
  print(utils::head(x$scores), ...)
  invisible(x)
}

#' @export
summary.al_result <- function(object, ...) {
  x <- object$scores$score
  structure(list(method = object$audit$method, scale = object$audit$scale,
    reference = object$reference$label, n = length(x), n_valid = sum(!is.na(x)),
    n_missing = sum(is.na(x)), score = summary(x),
    coverage = summary(object$scores$coverage),
    status = table(object$scores$status)), class = "summary.al_result")
}

#' @export
print.summary.al_result <- function(x, ...) {
  cat("Allostatic load summary\n")
  cat("Method:", x$method, "| scale:", x$scale, "| reference:", x$reference, "\n")
  cat("Records:", x$n, "| scored:", x$n_valid, "| not scored:", x$n_missing, "\n\nScore:\n")
  print(x$score, ...)
  cat("\nCoverage:\n"); print(x$coverage, ...)
  cat("\nStatus:\n"); print(x$status, ...)
  invisible(x)
}

#' @export
print.al_reference <- function(x, ...) {
  spec <- x$spec
  cat("Allostatic load reference:", x$label, "\n")
  cat("Method:", spec$method, "| biomarkers:", length(spec$biomarkers), "| balance:", spec$balance,
    if (length(spec$strata)) paste("| strata:", paste(spec$strata, collapse = ", ")),
    if (length(spec$age)) paste("| age model:", spec$age, "(df =", spec$age_df, ")"), "\n")
  print(utils::head(x$cutoffs, 12L), ...)
  if (nrow(x$cutoffs) > 12L) cat("...", nrow(x$cutoffs) - 12L, "more rows\n")
  invisible(x)
}

#' @export
as.data.frame.al_result <- function(x, row.names = NULL, optional = FALSE,
                                    what = c("scores", "components", "systems"), ...) {
  what <- match.arg(what)
  as.data.frame(x[[what]], row.names = row.names, optional = optional, ...)
}

#' @export
as.data.frame.al_survey <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$estimates, row.names = row.names, optional = optional, ...)
}

#' @export
as.data.frame.al_decomposition <- function(x, row.names = NULL, optional = FALSE,
                                           what = c("contributions", "summary"), ...) {
  what <- match.arg(what)
  as.data.frame(x[[what]], row.names = row.names, optional = optional, ...)
}

#' @export
as.data.frame.al_patterns <- function(x, row.names = NULL, optional = FALSE,
                                      what = c("patterns", "cooccurrence", "pattern_groups"), ...) {
  what <- match.arg(what)
  tbl <- x[[what]]
  if (is.null(tbl)) {
    stop("Component '", what, "' is not available in this al_patterns object.", call. = FALSE)
  }
  as.data.frame(tbl, row.names = row.names, optional = optional, ...)
}

#' @export
as.data.frame.al_bootstrap <- function(x, row.names = NULL, optional = FALSE,
                                       what = c("summary", "cutoff_intervals", "draws"), ...) {
  what <- match.arg(what)
  as.data.frame(x[[what]], row.names = row.names, optional = optional, ...)
}

#' @export
as.data.frame.al_multiverse <- function(x, row.names = NULL, optional = FALSE,
                                        what = c("summary", "effects", "stability", "agreement", "choices"), ...) {
  what <- match.arg(what)
  tbl <- x[[what]]
  if (is.null(tbl)) {
    stop("Component '", what, "' is not available in this al_multiverse object.", call. = FALSE)
  }
  as.data.frame(tbl, row.names = row.names, optional = optional, ...)
}

#' Simulate a small synthetic biomarker panel
#' @param n Number of participants.
#' @param seed Local random seed for reproducible results.
#' @param covariates When TRUE, adds `age` (20-85 years) and `sex`, with blood
#'   pressure, HDL and CRP depending on them, for normative-reference examples.
#'   The default output is unchanged from earlier versions.
#' @return A data frame with id, group, sbp, hdl and crp (plus age and sex).
#'   Values are synthetic and do not represent a real population.
#' @family result access
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
simulate_biomarkers <- function(n = 100L, seed = 1L, covariates = FALSE) {
  .check_integer(n, "n", 1L, text = "n must be a positive integer.")
  .check_seed(seed)
  .check_flag(covariates, "covariates")
  .with_seed(seed, function() {
    if (!covariates) return(data.frame(id = seq_len(n), group = rep(c("A", "B"), length.out = n),
      sbp = round(stats::rnorm(n, 125, 15), 1),
      hdl = round(pmax(10, stats::rnorm(n, 52, 12)), 1),
      crp = round(stats::rlnorm(n, 0, 0.8), 2)))
    age <- round(stats::runif(n, 20, 85), 1)
    sex <- sample(c("Female", "Male"), n, replace = TRUE)
    male <- sex == "Male"
    data.frame(id = seq_len(n), group = rep(c("A", "B"), length.out = n), age = age, sex = sex,
      sbp = round(stats::rnorm(n, ifelse(male, 113 + 0.42 * (age - 20), 104 + 0.62 * (age - 20)), 9 + 0.08 * (age - 20)), 1),
      hdl = round(pmax(10, stats::rnorm(n, ifelse(male, 47, 58) - 0.05 * (age - 20), 11)), 1),
      crp = round(stats::rlnorm(n, -0.4 + 0.012 * (age - 20), 0.8), 2))
  })
}
