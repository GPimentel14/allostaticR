#' Publication-oriented allostatic load theme
#' @param base_size Base font size in points.
#' @param base_family Font family available on the rendering device.
#' @return A ggplot2 theme without global side effects.
#' @family graphics
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
theme_allostatic <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1), hjust = 0, margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(face = "bold", colour = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(face = "bold", colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "black", size = ggplot2::rel(1)),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(colour = "black"),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(10, 12, 10, 10))
}

.plot_labels <- function(locale) {
  list(distribution = "Allostatic load distribution", biomarkers = "Biomarker risk profile",
    coverage = "Biomarker completeness", groups = "Score frequencies by group",
    boxplot = "Allostatic load by group", score = "Allostatic load", count = "Participants",
    risk = "At-risk participants (%)", observed = "Observed measurements (%)", group = "Group",
    missing = "Missing group", reference = "Reference", valid = "scored", absent = "not scored")
}

.score_axis <- function(p, scale, n_markers, axis = c("x", "y")) {
  axis <- match.arg(axis)
  continuous <- if (axis == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
  if (scale == "count") return(p + continuous(breaks = 0:n_markers, limits = c(0, n_markers),
    expand = ggplot2::expansion(mult = 0.04)))
  if (scale %in% c("proportion", "system_proportion")) return(p + continuous(limits = c(0, 1),
    labels = .percent, expand = ggplot2::expansion(mult = 0.04)))
  p
}

#' Generate automatic editable graphics
#' @param result An `al_result`.
#' @param locale Label language; only English (en) is supported.
#' @return Named list of ggplot objects with class al_plot_set. The omitted
#'   attribute records reasons for unavailable plots. Summaries use all rows.
#' @details Group plots require a group column retained during scoring. No files
#'   are written and no plots are printed. Biomarker denominators exclude NA;
#'   missing contributions are never counted as low risk. Count-score group
#'   frequencies use within-group proportions. Boxplots show medians,
#'   interquartile ranges and whiskers extending to 1.5 IQR, with all scored
#'   participants overlaid using reproducible categorical jitter. Without a group
#'   column, the boxplot summarizes the overall cohort. Cutoff panels use
#'   independent axes in the declared biomarker units; user labels are preserved.
#'   z-score protocols show mean z-scores per biomarker instead of risk
#'   prevalences. Age-conditional references omit the cutoff panel; use
#'   [plot_centiles()].
#' @family graphics
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @keywords internal
.plot_set <- function(result, locale = "en") {
  .assert_result(result); locale <- match.arg(locale, "en"); lab <- .plot_labels(locale)
  scores <- result$scores; cmp <- result$components; scale <- result$audit$scale
  count <- scale == "count"; zscale <- scale %in% c("z", "system_z")
  n_markers <- scores$n_total[1]
  score_label <- if (count) lab$score else .scale_label(scale)
  caption <- paste0(lab$reference, ": ", result$reference$label, " | ",
    result$audit$method, " | ", scale, " | n = ", nrow(scores),
    " | ", sum(!is.na(scores$score)), " ", lab$valid, " | ",
    sum(is.na(scores$score)), " ", lab$absent)
  wrap_labels <- function(x) vapply(x, function(y) paste(strwrap(y, 34), collapse = "\n"), character(1))
  plots <- list(); omitted <- list()
  valid <- scores[!is.na(scores$score), , drop = FALSE]
  if (nrow(valid)) {
    if (count) {
      freq <- data.frame(score = 0:n_markers)
      freq$n <- tabulate(match(valid$score, freq$score), nbins = nrow(freq))
      p <- ggplot2::ggplot(freq, ggplot2::aes(x = .data$score, y = .data$n)) +
        ggplot2::geom_col(fill = .al_palette$groups[["ink"]], width = 0.78) +
        ggplot2::scale_x_continuous(breaks = freq$score) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08)))
    } else if (!zscale) {
      p <- ggplot2::ggplot(valid, ggplot2::aes(x = .data$score)) +
        ggplot2::geom_histogram(binwidth = 0.1, boundary = 0, fill = .al_palette$groups[["ink"]], color = "white") +
        ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25), labels = .percent) +
        ggplot2::coord_cartesian(xlim = c(0, 1))
    } else {
      p <- ggplot2::ggplot(valid, ggplot2::aes(x = .data$score)) +
        ggplot2::geom_histogram(bins = 30, fill = .al_palette$groups[["ink"]], color = "white") +
        ggplot2::geom_vline(xintercept = 0, color = .al_palette$secondary, linewidth = 0.4)
    }
    plots$distribution <- .finish_plot(p, title = lab$distribution, caption = caption,
      x = score_label, y = lab$count,
      subtitle = sprintf("Scored participants: %s | Median: %.2f | Mean: %.2f",
        nrow(valid), stats::median(valid$score), mean(valid$score)))
  } else omitted$distribution <- "No valid scores."
  profile <- do.call(rbind, lapply(unique(cmp$biomarker), function(b) {
    d <- cmp[cmp$biomarker == b, ]; n <- sum(!is.na(d$contribution))
    data.frame(biomarker = b, label = paste0(d$label[1], " [", b, "]"),
      n_observed = n, proportion = if (n && !zscale) .mean_or_na(d$risk) else NA_real_,
      mean_z = if (n && zscale) mean(d$contribution, na.rm = TRUE) else NA_real_,
      observed = n / nrow(scores))
  }))
  profile$display <- paste0(profile$label, " (n=", profile$n_observed, ")")
  profile$display <- factor(profile$display, levels = rev(profile$display))
  if (!zscale && nrow(available <- profile[!is.na(profile$proportion), , drop = FALSE])) plots$biomarkers <- .finish_plot(
    ggplot2::ggplot(available, ggplot2::aes(x = .data$proportion, y = .data$display)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$proportion, yend = .data$display),
        color = .al_palette$groups[["light"]], linewidth = 1.4) +
      ggplot2::geom_point(color = .al_palette$groups[["ink"]], size = 3) +
      ggplot2::scale_x_continuous(limits = c(0, 1), labels = .percent) +
      ggplot2::scale_y_discrete(labels = wrap_labels),
    title = lab$biomarkers, caption = caption,
    x = lab$risk, y = NULL, subtitle = "Each denominator includes observed measurements only")
  else if (zscale && nrow(available <- profile[!is.na(profile$mean_z), , drop = FALSE])) plots$biomarkers <- .finish_plot(
    ggplot2::ggplot(available, ggplot2::aes(x = .data$mean_z, y = .data$display)) +
      ggplot2::geom_vline(xintercept = 0, color = .al_palette$muted, linewidth = 0.4) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$mean_z, yend = .data$display),
        color = .al_palette$groups[["light"]], linewidth = 1.4) +
      ggplot2::geom_point(color = .al_palette$groups[["ink"]], size = 3) +
      ggplot2::scale_x_continuous(limits = c(-1, 1) * max(0.5, 1.1 * max(abs(available$mean_z))),
        labels = function(x) formatC(x, format = "f", digits = 2)) +
      ggplot2::scale_y_discrete(labels = wrap_labels),
    title = lab$biomarkers, caption = caption,
    x = "Mean risk-oriented z-score (positive = towards risk)", y = NULL,
    subtitle = "Observed measurements only; reference mean = 0")
  else omitted$biomarkers <- "No observed biomarker values."
  plots$coverage <- .finish_plot(
    ggplot2::ggplot(profile, ggplot2::aes(x = .data$observed, y = .data$display)) +
      ggplot2::geom_col(fill = .al_palette$mid_grey, width = 0.65) +
      ggplot2::scale_x_continuous(limits = c(0, 1), labels = .percent,
        expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::scale_y_discrete(labels = wrap_labels),
    title = lab$coverage, caption = caption,
    x = lab$observed, y = NULL, subtitle = "Observed measurements as a percentage of all participants")
  if ("group" %in% names(valid) && nrow(valid)) {
    # Missing groups get their own category; avoid collision with actual labels.
    g <- as.character(valid$group); missing_label <- lab$missing
    while (missing_label %in% g[!is.na(g)]) missing_label <- paste0(missing_label, "*")
    g[is.na(g)] <- missing_label; valid$group <- g
    if (count) {
      tab <- as.data.frame(table(group = g, score = factor(valid$score, levels = 0:n_markers)), stringsAsFactors = FALSE)
      tab$score <- as.numeric(as.character(tab$score))
      tab$proportion <- tab$Freq / stats::ave(tab$Freq, tab$group, FUN = sum)
      p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$score, y = .data$proportion)) +
        ggplot2::geom_col(fill = .al_palette$groups[["ink"]], width = 0.78) +
        ggplot2::facet_wrap(~group, labeller = ggplot2::label_wrap_gen(width = 26)) +
        ggplot2::scale_y_continuous(labels = .percent) +
        ggplot2::scale_x_continuous(breaks = 0:n_markers)
    } else {
      p <- ggplot2::ggplot(valid, ggplot2::aes(x = .data$group, y = .data$score)) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 105),
          colour = "grey45", alpha = 0.35, size = 1) +
        ggplot2::geom_boxplot(fill = NA, colour = "black", width = 0.5, outlier.shape = NA, linewidth = 0.5)
      if (!zscale) p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), labels = .percent)
    }
    group_n <- table(valid$group)
    plots$groups <- .finish_plot(p, title = if (count) lab$groups else lab$boxplot, caption = caption,
      x = if (count) lab$score else lab$group, y = if (count) "Proportion of scored group" else score_label,
      subtitle = paste(paste0(names(group_n), ": n=", as.integer(group_n)), collapse = " | "))
  } else omitted$groups <- "No grouping column or no valid scores."
  if (nrow(valid)) {
    boxes <- valid
    if (!"group" %in% names(boxes)) boxes$group <- "All participants"
    g <- as.character(boxes$group)
    missing_label <- lab$missing
    while (missing_label %in% g[!is.na(g)]) missing_label <- paste0(missing_label, "*")
    g[is.na(g)] <- missing_label
    boxes$group <- factor(g, levels = rev(unique(g)))
    group_n <- table(boxes$group)
    p <- ggplot2::ggplot(boxes, ggplot2::aes(x = .data$score, y = .data$group)) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0, height = 0.15, seed = 104),
        colour = "grey45", alpha = 0.35, size = 1) +
      ggplot2::geom_boxplot(orientation = "y", fill = NA, colour = "black", width = 0.5,
        outlier.shape = NA, linewidth = 0.5) +
      ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(y)
        paste(strwrap(paste0(y, " (n=", group_n[y], ")"), 30), collapse = "\n"), character(1)))
    plots$boxplot <- .finish_plot(.score_axis(p, scale, n_markers, "x"),
      title = if ("group" %in% names(scores)) lab$boxplot else "Allostatic load summary", caption = caption,
      x = score_label, y = NULL,
      subtitle = "Box: median and interquartile range | diamond: mean | points: scored participants")
  } else omitted$boxplot <- "No valid scores."
  if (length(result$reference$spec$age)) omitted$cutoffs <- "Thresholds vary with age; use plot_centiles()."
  else if (any(is.finite(cmp$value))) plots$cutoffs <- plot_cutoffs(result)
  else omitted$cutoffs <- "No observed biomarker values."
  structure(plots, class = c("al_plot_set", "list"), omitted = omitted)
}

#' Select an allostatic load graphic
#' @param result An `al_result`.
#' @param type Graphic name: distribution, biomarkers, coverage, groups,
#'   boxplot or cutoffs.
#' @param locale Label language; only English (en) is supported.
#' @return An editable ggplot object. Customize with ggplot2 labs(), scales,
#'   coordinates and themes. To choose all geometries and mappings yourself,
#'   use plot_data() with ggplot2::ggplot().
#' @family graphics
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @export
plot_load <- function(result, type = c("distribution", "biomarkers", "coverage", "groups", "boxplot", "cutoffs"),
                     locale = "en") {
  type <- match.arg(type); plots <- .plot_set(result, match.arg(locale, "en"))
  if (is.null(plots[[type]])) stop("Graphic unavailable: ", attr(plots, "omitted")[[type]], call. = FALSE)
  plots[[type]]
}

#' @export
plot.al_result <- function(x, ...) {
  p <- plot_load(x, ...); print(p); invisible(p)
}

#' Visualize measured biomarkers with fitted cuts
#' @param result An `al_result` retaining measured components.
#' @return Editable faceted histogram with threshold lines in declared units:
#'   risk thresholds for count methods, or the reference mean (geometric mean
#'   for log-transformed biomarkers) for z_mean. Stratified references draw one
#'   line per stratum. Medication-adjusted risks are not derivable from
#'   threshold lines alone. Age-conditional references use [plot_centiles()].
#' @family graphics
#' @seealso [allostaticR-package][] for the package overview and workflow.
#' @keywords internal
plot_cutoffs <- function(result) {
  .assert_result(result)
  spec <- result$reference$spec
  if (length(spec$age)) stop("Thresholds vary with age; use plot_centiles().", call. = FALSE)
  facet_of <- function(d) paste0(d$label, " [", d$biomarker, "] (", d$unit, ")")
  d <- result$components; d$facet <- facet_of(d)
  cuts <- result$reference$cutoffs; cuts$facet <- facet_of(cuts)
  stratum <- if (length(spec$strata)) cuts$stratum else rep("all", nrow(cuts))
  if (spec$method == "z_mean") {
    transform <- .biomarker_field(spec, "transform")[match(cuts$biomarker,
      .biomarker_field(spec, "name"))]
    lines <- data.frame(facet = cuts$facet, stratum = stratum,
      cutoff = ifelse(transform == "log", exp(cuts$center), cuts$center))
  } else lines <- rbind(data.frame(facet = cuts$facet, stratum = stratum, cutoff = cuts$lower),
                        data.frame(facet = cuts$facet, stratum = stratum, cutoff = cuts$upper))
  lines <- lines[is.finite(lines$cutoff), ]
  d <- d[is.finite(d$value), ]
  if (!nrow(d)) stop("No observed biomarker values to plot.", call. = FALSE)
  line_layer <- if (length(spec$strata))
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = .data$cutoff, linetype = .data$stratum),
      color = .al_palette$diverging[["high"]], inherit.aes = FALSE)
  else ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = .data$cutoff),
      color = .al_palette$diverging[["high"]], linetype = "dashed", inherit.aes = FALSE)
  what <- if (spec$method == "z_mean") "reference means" else "fixed risk thresholds"
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 20, fill = .al_palette$groups[["ink"]], color = "white") +
    ggplot2::scale_x_continuous(n.breaks = 4) + line_layer +
    ggplot2::facet_wrap(~facet, scales = "free", labeller = ggplot2::label_wrap_gen(width = 28))
  args <- list(title = if (spec$method == "z_mean") "Biomarker distributions and reference means"
      else "Biomarker distributions and reference thresholds",
    subtitle = paste("Independent panel axes | Lines mark", what), x = "Measured value (unit in panel title)",
    y = "Observed participants (count)",
    caption = paste("Reference:", result$reference$label,
      "| lines:", what, "; missing measurements excluded. Medication rules, if present, additionally affect risk."))
  if (length(spec$strata)) args$linetype <- "Stratum"
  do.call(.finish_plot, c(list(p), args))
}
