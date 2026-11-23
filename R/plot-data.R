#' Extract data for user-defined ggplot2 graphics
#' @param result An `al_result` from score_load() or allostatic_load().
#' @param type Data to return: scores (one row per participant), biomarkers
#'   (one row per participant and biomarker), systems (one row per participant
#'   and system), or cutoffs (one row per biomarker).
#' @return A plain data frame, without selecting a geometry, axes or theme.
#'   Missing observations are retained. Group metadata, when available, is
#'   attached to participant-level tables in the original scoring order.
#' @details Biomarker values retain their declared units; use separate panels
#'   and scales when plotting biomarkers with different units. System scores
#'   are observed risk proportions, not full-panel allostatic load counts.
#' @examples
#' source(system.file("examples", "panel.R", package = "allostaticR"))
#' result <- allostatic_load(data, spec, id = "id", group = "group")
#' scores <- plot_data(result)
#' # Count only participants with a valid score in the group labels.
#' group_n <- table(scores$group[!is.na(scores$score)])
#' group_labels <- paste0(names(group_n), "\nn = ", as.integer(group_n))
#' names(group_labels) <- names(group_n)
#'
#' p <- ggplot2::ggplot(scores, ggplot2::aes(x = group, y = score)) +
#'   ggplot2::geom_boxplot(
#'     width = 0.42, fill = "#E4EEF5", colour = "#24465C",
#'     linewidth = 0.65, outlier.shape = NA, na.rm = TRUE
#'   ) +
#'   ggplot2::geom_point(
#'     position = ggplot2::position_jitter(width = 0.09, height = 0, seed = 42),
#'     colour = "#27658A", size = 1.9, alpha = 0.55, na.rm = TRUE
#'   ) +
#'   ggplot2::scale_x_discrete(labels = group_labels) +
#'   ggplot2::scale_y_continuous(
#'     breaks = 0:2, limits = c(0, 2),
#'     expand = ggplot2::expansion(add = 0.18)
#'   ) +
#'   ggplot2::labs(
#'     title = "Allostatic load by study group",
#'     subtitle = "Two-biomarker panel: systolic blood pressure and HDL cholesterol",
#'     x = "Study group", y = "At-risk biomarkers (count)",
#'     caption = "Synthetic data. Boxes: median and IQR. Points: individual participants."
#'   ) +
#'   ggplot2::theme_classic(base_size = 12) +
#'   ggplot2::theme(
#'     plot.title = ggplot2::element_text(face = "bold", size = 15, colour = "#213547"),
#'     plot.subtitle = ggplot2::element_text(size = 10, colour = "#526574",
#'       margin = ggplot2::margin(b = 16)),
#'     plot.caption = ggplot2::element_text(size = 9, colour = "#526574", hjust = 0,
#'       margin = ggplot2::margin(t = 12)),
#'     plot.title.position = "plot", plot.caption.position = "plot",
#'     axis.text = ggplot2::element_text(colour = "#334155"),
#'     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
#'     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
#'     axis.line = ggplot2::element_line(colour = "#94A3B8", linewidth = 0.4),
#'     axis.ticks = ggplot2::element_line(colour = "#94A3B8"),
#'     plot.margin = ggplot2::margin(14, 16, 12, 12)
#'   )
#' p
#' @family graphics
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
plot_data <- function(result, type = c("scores", "biomarkers", "systems", "cutoffs")) {
  .assert_result(result)
  type <- match.arg(type)
  data <- switch(type, scores = result$scores, biomarkers = result$components,
    systems = result$systems, cutoffs = result$reference$cutoffs)
  if (type %in% c("biomarkers", "systems") && "group" %in% names(result$scores))
    data$group <- result$scores$group[match(data$row, result$scores$row)]
  data
}
