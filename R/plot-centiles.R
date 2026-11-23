#' Plot age-conditional reference curves
#'
#' Draws growth-chart style curves of an age-conditional reference: the
#' age-specific risk thresholds (count methods) or the reference mean and
#' mean +/- 1.96 SD (z_mean) for each biomarker and stratum, in declared units.
#' With a scored result, measured values are overlaid and coloured by risk.
#'
#' @param x An `al_result` scored against an age-conditional reference, or the
#'   `al_reference` itself (curves only).
#' @param biomarkers Optional biomarker names to show.
#' @return An editable ggplot object. Curves span the observed reference age
#'   range of each stratum; points outside it are flagged, never extrapolated.
#' @family graphics
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
plot_centiles <- function(x, biomarkers = NULL) {
  reference <- if (inherits(x, "al_result")) x$reference else x
  .assert_reference(reference); spec <- reference$spec
  if (!length(spec$age)) stop("plot_centiles() requires an age-conditional reference.", call. = FALSE)
  all_names <- .biomarker_field(spec, "name")
  if (is.null(biomarkers)) biomarkers <- all_names
  if (!is.character(biomarkers) || !length(biomarkers) || anyNA(biomarkers) || !all(biomarkers %in% all_names))
    stop("biomarkers must name protocol biomarkers.", call. = FALSE)
  facet_label <- function(name) {
    b <- spec$biomarkers[match(name, all_names)]
    factor(paste0(vapply(b, `[[`, character(1), "label"), " (", vapply(b, `[[`, character(1), "unit"), ")"),
      levels = unique(paste0(vapply(spec$biomarkers[match(biomarkers, all_names)], `[[`, character(1), "label"), " (",
        vapply(spec$biomarkers[match(biomarkers, all_names)], `[[`, character(1), "unit"), ")")))
  }
  models <- Filter(function(m) m$biomarker %in% biomarkers, reference$models)
  curves <- do.call(rbind, lapply(models, function(m) {
    ages <- seq(m$age_range[1], m$age_range[2], length.out = 120L); G <- .spline_design(ages, m)
    if (spec$method == "z_mean") {
      center <- drop(G %*% m$coef$center); scale <- drop(G %*% m$coef$scale)
      do.call(rbind, lapply(c(-1.96, 0, 1.96), function(z) data.frame(biomarker = m$biomarker, stratum = m$stratum,
        age = ages, value = .untransform(center + z * scale, m$transform), z = z,
        curve = if (z == 0) "Reference mean" else "Mean \u00b1 1.96 SD")))
    } else do.call(rbind, lapply(names(m$coef), function(p) data.frame(biomarker = m$biomarker, stratum = m$stratum,
        age = ages, value = .untransform(drop(G %*% m$coef[[p]]), m$transform), z = if (p == "lower") -1 else 1,
        curve = paste(if (p == "lower") "Lower" else "Upper", "risk threshold"))))
  }))
  curves$facet <- facet_label(curves$biomarker)
  curves$path <- paste(curves$stratum, curves$curve, curves$z)
  p <- ggplot2::ggplot()
  if (inherits(x, "al_result")) {
    cmp <- x$components[x$components$biomarker %in% biomarkers & is.finite(x$components$value), , drop = FALSE]
    cmp$age <- x$scores$age[cmp$row]
    if (!length(spec$strata)) cmp$stratum <- "all"
    cmp$facet <- facet_label(cmp$biomarker)
    cmp$class <- if (.is_count_method(spec)) ifelse(is.na(cmp$risk), "Outside reference age",
      ifelse(cmp$risk == 1, paste("At risk:", cmp$system), "Not at risk")) else
      ifelse(is.na(cmp$contribution), "Outside reference age", "Measured")
    colors <- .system_colors(.protocol_systems(spec))
    values <- c(stats::setNames(colors, paste("At risk:", names(colors))), `Not at risk` = .al_palette$light_grey,
      Measured = .al_palette$grey, `Outside reference age` = .al_palette$diverging[["high"]])
    p <- p + ggplot2::geom_point(data = cmp, ggplot2::aes(x = .data$age, y = .data$value, colour = .data$class),
        size = 0.9, alpha = 0.65) +
      ggplot2::scale_colour_manual(values = values[names(values) %in% unique(cmp$class)])
  }
  p <- p + ggplot2::geom_line(data = curves, ggplot2::aes(x = .data$age, y = .data$value, group = .data$path,
      linetype = .data$curve), colour = .al_palette$primary, linewidth = 0.7) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22", "solid", "42"),
      c("Reference mean", "Mean \u00b1 1.96 SD", "Upper risk threshold", "Lower risk threshold")))
  p <- p + if (length(spec$strata)) ggplot2::facet_grid(facet ~ stratum, scales = "free_y",
      labeller = ggplot2::label_wrap_gen(width = 14))
    else ggplot2::facet_wrap(~facet, scales = "free_y", labeller = ggplot2::label_wrap_gen(width = 28))
  p <- .finish_plot(p,
    title = "Age-conditional reference curves",
    subtitle = if (.is_count_method(spec)) sprintf("Age-specific risk thresholds (quantile regression, natural spline df = %s)", spec$age_df)
      else sprintf("Age-specific mean and SD (absolute-residual method, natural spline df = %s)", spec$age_df),
    x = paste0("Age (", spec$age, ")"), y = "Measured value", colour = NULL, linetype = NULL,
    caption = paste0("Reference: ", reference$label, ". Curves span the observed reference ages",
      if (length(spec$strata)) " of each stratum" else "", "; values outside are not extrapolated",
      if (spec$age_extrapolate) " unless the protocol allows it" else "", ".")) +
    ggplot2::theme(legend.box = "vertical", legend.box.just = "left") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, order = 2, override.aes = list(size = 2.5, alpha = 1)),
      linetype = ggplot2::guide_legend(order = 1))
}
