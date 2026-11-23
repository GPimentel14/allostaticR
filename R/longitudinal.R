#' Score repeated visits with a fixed reference
#' @param data Long-format biomarker data frame.
#' @param reference Fixed reference fitted independently or at baseline.
#' @param id,time Participant and visit-time columns. Times must be numeric,
#'   Date or POSIXct, finite and non-missing. Duplicate participant/time pairs
#'   are rejected. Each row remains a distinct scored visit.
#' @param time_unit Label for numeric time units; dates use days in visit_changes().
#' @param subjects Optional identifiers to show individually (at most 20).
#'   By default only cohort mean trajectories are plotted; their composition
#'   can change across visits and they are not paired change estimates.
#' @param group Optional column with a constant group per participant; cohort
#'   means are then summarized and plotted by group.
#' @return An `al_trajectory` containing visit scores, visit components,
#'   per-time coverage and means, a fixed reference and editable plots. No
#'   mixed model or causal effect is fitted.
#' @family longitudinal analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
score_trajectory <- function(data, reference, id, time, time_unit = "visit", subjects = NULL, group = NULL) {
  .assert_reference(reference); validate_data(data, reference$spec)
  .scalar_string(id, "id"); .scalar_string(time, "time"); .scalar_string(time_unit, "time_unit")
  if (!all(c(id, time) %in% names(data))) stop("Missing id/time column.", call. = FALSE)
  ids <- data[[id]]; times <- data[[time]]
  if (!is.atomic(ids) || !is.null(dim(ids)) || anyNA(ids)) stop("Participant ids must be non-missing vectors.", call. = FALSE)
  if ((!is.numeric(times) && !inherits(times, c("Date", "POSIXct"))) ||
      !is.null(dim(times)) || any(!is.finite(as.numeric(times))))
    stop("Visit times must be finite numeric, Date or POSIXct values.", call. = FALSE)
  if (anyDuplicated(data.frame(id = ids, time = times))) stop("Duplicate participant/time pairs.", call. = FALSE)
  if (!is.null(group)) {
    .scalar_string(group, "group")
    if (!group %in% names(data)) stop("Missing group column: ", group, call. = FALSE)
    if (all(is.na(data[[group]]))) stop("Group column has no non-missing values.", call. = FALSE)
    per_person <- tapply(as.character(data[[group]]), ids, function(g) length(unique(g)))
    if (any(per_person > 1L)) stop("Each participant must keep one group across visits.", call. = FALSE)
  }
  result <- score_load(data, reference, group = group)
  visits <- result$scores; visits$id <- ids; visits$time <- times
  components <- result$components; components$id <- ids[components$row]; components$time <- times[components$row]
  visits <- visits[order(visits$id, visits$time), , drop = FALSE]
  unique_times <- sort(unique(times))
  grouped <- !is.null(group)
  keys <- if (grouped) unique(stats::na.omit(as.character(visits$group))) else "All participants"
  summary <- do.call(rbind, lapply(keys, function(g) do.call(rbind, lapply(seq_along(unique_times), function(i) {
    d <- visits[visits$time == unique_times[i] & (if (grouped) as.character(visits$group) %in% g else TRUE), ]
    out <- data.frame(time = unique_times[i], n = nrow(d), n_valid = sum(is.finite(d$score)),
      mean_score = .mean_or_na(d$score))
    if (grouped) out$group <- g
    out
  }))))
  caption <- paste0("Fixed reference: ", reference$label, "; scale: ", result$audit$scale,
    ". Cohort means use available visits; composition may vary.")
  aes_line <- if (grouped) ggplot2::aes(x = .data$time, y = .data$mean_score, colour = .data$group, group = .data$group)
    else ggplot2::aes(x = .data$time, y = .data$mean_score)
  line_args <- list(linewidth = 0.8, na.rm = TRUE); point_args <- list(size = 2.5, na.rm = TRUE)
  if (!grouped) line_args$colour <- point_args$colour <- .al_palette$groups[["ink"]]
  cohort_args <- list(title = "Allostatic load over time",
    subtitle = paste0("Participants: ", length(unique(ids)), " | Valid visits per time: ",
      min(summary$n_valid), " to ", max(summary$n_valid)), x = if (inherits(times, c("Date", "POSIXct"))) "Visit date" else paste("Time (", time_unit, ")", sep = ""),
    y = paste("Mean score (", result$audit$scale, ")", sep = ""), caption = caption)
  if (grouped) cohort_args$colour <- "Group"
  p <- do.call(.finish_plot, c(list(ggplot2::ggplot(summary, aes_line) +
    (if (length(unique_times) > 1L) do.call(ggplot2::geom_line, line_args) else NULL) +
    do.call(ggplot2::geom_point, point_args)), cohort_args))
  if (grouped) p <- p +
    ggplot2::scale_colour_manual(values = .palette_for(length(keys), .al_palette$categorical))
  plots <- list(cohort = p)
  if (!is.null(subjects)) {
    if (!length(subjects) || anyNA(subjects) || anyDuplicated(subjects) || length(subjects) > 20 || !all(subjects %in% ids))
      stop("subjects must select 1 to 20 unique existing participants.", call. = FALSE)
    d <- visits[visits$id %in% subjects, ]; d$participant <- as.character(d$id)
    plots$individual <- .finish_plot(
      ggplot2::ggplot(d,
          ggplot2::aes(x = .data$time, y = .data$score, group = .data$participant)) +
        ggplot2::geom_line(color = .al_palette$groups[["ink"]], na.rm = TRUE) + ggplot2::geom_point(color = .al_palette$groups[["ink"]], na.rm = TRUE) +
        ggplot2::facet_wrap(~participant),
      title = "Selected participant trajectories", x = "Visit time", y = paste("Score (", result$audit$scale, ")", sep = ""),
      caption = "Source: supplied visits. Same reference for every visit; missing scores remain missing.")
  }
  structure(list(visits = visits, components = components, summary = summary, reference = reference,
    time_unit = if (inherits(times, c("Date", "POSIXct"))) "days" else time_unit,
    scale = result$audit$scale, plots = plots), class = "al_trajectory")
}

#' Calculate paired first-to-last visit changes
#' @param trajectory Object from score_trajectory().
#' @param variation Optional biological-variation table (see
#'   [measurement_uncertainty()]). For count protocols, each first-to-last
#'   change in a biomarker's risk flag is then classified with the reference
#'   change value (RCV): a flip is `confirmed` when the measured change exceeds
#'   what analytical and within-person biological variation explain.
#' @param conf Two-sided confidence level of the RCV.
#' @return One row per participant with first and last observed visit, scores,
#'   coverage, elapsed time, change and reason. A missing endpoint remains
#'   missing; intermediate valid visits never silently replace endpoints. With
#'   `variation`, columns `n_flips`, `n_confirmed_flips` and `confirmed_change`
#'   are added and attribute `biomarker_changes` lists every biomarker change.
#' @references Fokkema MR et al. (2006). Lognormal reference change values.
#'   Clinical Chemistry 52, 1602-1603.
#' @family longitudinal analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
visit_changes <- function(trajectory, variation = NULL, conf = 0.95) {
  if (!inherits(trajectory, "al_trajectory")) stop("Expected an score_trajectory object.", call. = FALSE)
  visits <- trajectory$visits
  split_visits <- split(seq_len(nrow(visits)), factor(visits$id, levels = unique(visits$id)))
  rows <- lapply(split_visits, function(i) {
    d <- visits[i, ]; first <- d[1, ]; last <- d[nrow(d), ]
    elapsed <- if (inherits(d$time, c("Date", "POSIXct")))
      as.numeric(difftime(last$time, first$time, units = "days")) else last$time - first$time
    status <- if (nrow(d) < 2L) "single_visit" else if (anyNA(c(first$score, last$score))) "missing_endpoint" else "ok"
    data.frame(id = first$id, first_time = first$time, last_time = last$time,
      first_score = first$score, last_score = last$score,
      first_coverage = first$coverage, last_coverage = last$coverage,
      elapsed = elapsed, time_unit = trajectory$time_unit,
      change = if (status == "ok") last$score - first$score else NA_real_, status = status)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (is.null(variation)) return(out)
  .rcv_changes(trajectory, out, variation, conf)
}

#' Model allostatic load trajectories with a linear mixed model
#'
#' Uses every valid visit of a trajectory, not only its endpoints: the score is
#' modelled as a linear function of time (by group when the trajectory was
#' scored with a group column) with participant random intercepts, and
#' optionally random slopes, fitted by restricted maximum likelihood with
#' nlme. Count scores are treated as approximately continuous.
#'
#' @param trajectory Object from score_trajectory().
#' @param random `"slope"` (random intercept and slope) or `"intercept"`.
#' @param conf Confidence level for fixed effects and mean curves.
#' @return An `al_trajectory_model` with `fixed` (estimates, standard errors
#'   and intervals; time in the trajectory's units, years for dates),
#'   `random` (participant-level SDs and residual SD), `random_structure`
#'   (`"slope"` or `"intercept"` actually fitted), `curves` (model-based
#'   mean trajectories with intervals), the nlme `model` and an editable plot
#'   overlaying observed visit means. Associations are descriptive; missing
#'   visits are assumed missing at random.
#' @family longitudinal analysis
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
model_trajectory <- function(trajectory, random = c("slope", "intercept"), conf = 0.95) {
  if (!inherits(trajectory, "al_trajectory")) stop("Expected an score_trajectory object.", call. = FALSE)
  random <- match.arg(random); .check_open_unit(conf, "conf")
  .require_suggested("nlme", "fit mixed-effects trajectory models")
  v <- trajectory$visits[!is.na(trajectory$visits$score), , drop = FALSE]
  dated <- inherits(v$time, c("Date", "POSIXct"))
  v$t <- if (dated) as.numeric(difftime(v$time, min(v$time), units = "days")) / 365.25 else as.numeric(v$time)
  unit <- if (dated) "years since first visit" else trajectory$time_unit
  grouped <- "group" %in% names(v) && length(unique(stats::na.omit(v$group))) > 1L
  if (grouped) { v <- v[!is.na(v$group), ]; v$group <- factor(v$group) }
  if (length(unique(v$id)) < 3L || length(unique(v$t)) < 2L) stop("Need at least three participants and two distinct times.", call. = FALSE)
  v$id <- factor(v$id)
  fixed <- if (grouped) score ~ t * group else score ~ t
  re <- if (random == "slope") ~ t | id else ~ 1 | id
  actual_random <- random
  model <- tryCatch(nlme::lme(fixed, random = re, data = v, method = "REML",
      control = nlme::lmeControl(opt = "optim", returnObject = TRUE)),
    error = function(e) {
      if (random == "slope") {
        actual_random <<- "intercept"
        warning("Random-slope model failed; falling back to random intercepts only: ", conditionMessage(e), call. = FALSE)
        nlme::lme(fixed, random = ~ 1 | id, data = v, method = "REML")
      } else stop(e)
    })
  z <- stats::qnorm(1 - (1 - conf) / 2)
  tab <- summary(model)$tTable
  fixed_table <- data.frame(term = gsub("(^|:)t($|:)", "\\1time\\2", rownames(tab)), estimate = tab[, "Value"],
    std_error = tab[, "Std.Error"], lower = tab[, "Value"] - z * tab[, "Std.Error"], upper = tab[, "Value"] + z * tab[, "Std.Error"],
    df = tab[, "DF"])
  rownames(fixed_table) <- NULL
  vc <- nlme::VarCorr(model)
  random_table <- data.frame(component = rownames(vc), sd = suppressWarnings(as.numeric(vc[, "StdDev"])))
  grid <- expand.grid(t = seq(min(v$t), max(v$t), length.out = 50), group = if (grouped) levels(v$group) else "All participants")
  if (grouped) grid$group <- factor(grid$group, levels(v$group))
  X <- stats::model.matrix(stats::delete.response(stats::terms(fixed)), if (grouped) grid else grid["t"])
  beta <- nlme::fixef(model); V <- stats::vcov(model)
  grid$fit <- drop(X %*% beta); se <- sqrt(rowSums((X %*% V) * X))
  grid$lower <- grid$fit - z * se; grid$upper <- grid$fit + z * se
  v$grp <- if (grouped) as.character(v$group) else "All participants"
  observed <- stats::aggregate(score ~ t + grp, data = v, FUN = mean)
  names(observed)[2] <- "group"
  colors <- .palette_for(length(unique(grid$group)), .al_palette$categorical)
  plot <- .finish_plot(
    ggplot2::ggplot(grid, ggplot2::aes(x = .data$t, y = .data$fit, colour = .data$group, fill = .data$group)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.15, colour = NA) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(data = observed, ggplot2::aes(x = .data$t, y = .data$score, colour = .data$group), inherit.aes = FALSE, size = 2) +
      ggplot2::scale_colour_manual(values = colors) + ggplot2::scale_fill_manual(values = colors, guide = "none"),
    title = "Modelled allostatic load trajectories",
    subtitle = sprintf("Linear mixed model, random %s; lines: mean trajectory with %g%% interval; points: observed means",
      if (actual_random == "slope") "intercepts and slopes" else "intercepts", 100 * conf),
    x = paste0("Time (", unit, ")"), y = paste("Score (", trajectory$scale, ")", sep = ""), colour = NULL, width = 100,
    caption = sprintf("%d participants, %d valid visits. Fixed reference: %s.", length(unique(v$id)), nrow(v),
      trajectory$reference$label))
  structure(list(fixed = fixed_table, random = random_table, random_structure = actual_random, curves = grid,
    model = model, time_unit = unit, plots = list(trajectory = plot)), class = "al_trajectory_model")
}

#' @export
print.al_trajectory_model <- function(x, ...) {
  cat("Mixed-effects trajectory model (time:", x$time_unit, ")\n\n")
  print(x$fixed, row.names = FALSE, digits = 3, ...)
  cat("\n"); print(x$random, row.names = FALSE, digits = 3)
  invisible(x)
}
