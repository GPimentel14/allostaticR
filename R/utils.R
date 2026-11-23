# Shared internal helpers.

.package_version <- function() as.character(utils::packageVersion("allostaticR"))

.with_seed <- function(seed, f) {
  .check_seed(seed)
  present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (present) old <- get(".Random.seed", .GlobalEnv)
  on.exit(if (present) assign(".Random.seed", old, .GlobalEnv)
          else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed); f()
}

# A local random seed: non-negative integer within the RNG range.
.check_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be a non-negative integer.", call. = FALSE)
}

.check_integer <- function(x, name, min = 1L, max = Inf, text = NULL) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) || x < min || x > max)
    stop(if (is.null(text)) paste0(name, " must be an integer from ", min,
      if (is.finite(max)) paste0(" to ", max) else " upwards", ".") else text, call. = FALSE)
}

.check_open_unit <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop(name, " must be in (0, 1).", call. = FALSE)
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(name, " must be TRUE or FALSE.", call. = FALSE)
}

# A non-empty data frame whose column names are unique.
.check_data_frame <- function(data, name = "data") {
  if (!is.data.frame(data) || !nrow(data)) stop(name, " must be a non-empty data frame.", call. = FALSE)
  if (anyDuplicated(names(data))) stop("Column names must be unique.", call. = FALSE)
}

# Medication columns that are absent, not logical or hold missing values.
# Returns the offending names; the caller decides whether to abort or report.
.invalid_medication_columns <- function(data, columns) {
  columns <- unname(columns)
  if (!length(columns)) return(character())
  valid <- vapply(columns, function(column)
    column %in% names(data) && is.logical(data[[column]]) && !anyNA(data[[column]]), logical(1))
  columns[!valid]
}

.require_suggested <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the '", package, "' package to ", purpose, ".", call. = FALSE)
}

.valid_names <- function(x) !is.null(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)

.wrap_caption <- function(text, width = 85) paste(strwrap(text, width), collapse = "\n")

.percent <- function(x) paste0(round(100 * x), "%")

# Validated categorical order (colour-vision checked); systems keep their colour in every figure.
.al_palette <- list(
  categorical = c("#1f4e79", "#9b2226", "#2a6f62", "#b07d1a", "#5b4b8a", "#3d5a2a", "#a0522d", "#4d4d4d"),
  primary = "#0b0b0b", secondary = "#52514e", muted = "#898781", grid = "#e1e0d9",
  neutral = "#e6e6e6", unavailable = "#5b5a56",
  groups = c(ink = "#262626", mid = "#7f7f7f", light = "#bfbfbf"),
  diverging = c(low = "#1f4e79", mid = "#f2f2f2", high = "#9b2226"),
  light_grey = "#a6a6a6", light_blue = "#d9e2ec", mid_grey = "#595959",
  grey = "#8c8c8c", light_red = "#efd9d9")

.system_colors <- function(systems) {
  systems <- unique(as.character(systems))
  cols <- c(.al_palette$categorical, rep(.al_palette$muted, max(0L, length(systems) - 8L)))
  stats::setNames(cols[seq_along(systems)], systems)
}

.protocol_systems <- function(spec) unique(.biomarker_field(spec, "system"))

# One metadata field of every protocol biomarker, in protocol order.
.biomarker_field <- function(spec, field) vapply(spec$biomarkers, `[[`, character(1), field)

# Protocol metadata of every biomarker, in protocol order.
.biomarker_info <- function(spec) {
  data.frame(name = .biomarker_field(spec, "name"),
    label = .biomarker_field(spec, "label"),
    system = .biomarker_field(spec, "system"))
}

# One row per protocol biomarker, built by `build`; the helper adds no columns.
.per_biomarker <- function(spec, build) {
  out <- do.call(rbind, lapply(spec$biomarkers, build))
  rownames(out) <- NULL
  out
}

.is_count_method <- function(spec) spec$method != "z_mean"

.score_scale <- function(spec) {
  if (spec$method == "z_mean") return(if (spec$balance == "system") "system_z" else "z")
  if (spec$balance == "system") return("system_proportion")
  if (spec$missing == "complete") "count" else "proportion"
}

.scale_label <- function(scale) {
  switch(scale, count = "Allostatic load (count)", proportion = "At-risk proportion",
    system_proportion = "System-balanced risk proportion", z = "Mean risk-oriented z-score",
    system_z = "System-balanced mean z-score", scale)
}

# Risk fraction on a 0-1 scale; z scales are not risk fractions.
.risk_fraction <- function(result) {
  x <- result$scores$score
  switch(result$audit$scale, count = x / result$scores$n_total[1L], proportion = x,
    system_proportion = x, rep(NA_real_, length(x)))
}

# Comparable score: risk fraction for count scales, native units otherwise.
.comparable_score <- function(result) {
  if (result$audit$scale == "count") result$scores$score / result$scores$n_total[1L] else result$scores$score
}

.mean_or_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_

.safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(NA_real_)
  unname(stats::cor(x[ok], y[ok], method = "spearman"))
}

.percentile_interval <- function(x, conf) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  unname(stats::quantile(x, c((1 - conf) / 2, 1 - (1 - conf) / 2)))
}

# Single rule for risk flags, shared by scoring and by the imputation check.
.risk_flag <- function(values, p, direction, inclusive) {
  flag <- rep(FALSE, length(values))
  if (direction %in% c("low", "both")) flag <- flag | if (inclusive) values <= p$lower else values < p$lower
  if (direction %in% c("high", "both")) flag <- flag | if (inclusive) values >= p$upper else values > p$upper
  flag
}

# Mid-rank percentile; tied values share one rank.
.midrank_percentile <- function(x) {
  if (!length(x)) return(NA_real_)
  (rank(x) - 0.5) / length(x)
}

# Paired comparison of two score vectors.
.paired_scores <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  list(n_matched = sum(ok),
    mean_difference = .mean_or_na(y[ok] - x[ok]),
    mean_absolute_difference = .mean_or_na(abs(y[ok] - x[ok])),
    spearman = .safe_spearman(x, y))
}

# Colours for n elements, recycled from the given palette.
.palette_for <- function(n, colours = .al_palette$categorical) rep(colours, length.out = n)

# A reproducible random subset; the caller's random state is restored.
.random_subset <- function(x, size, seed = 1L) {
  if (size >= length(x)) return(x)
  x[.with_seed(seed, function() sort(sample.int(length(x), size)))]
}

# Print a table without row names and with a fixed number of digits.
.print_table <- function(x, digits = 3, ...) print(x, row.names = FALSE, digits = digits, ...)

# Median, minimum and maximum of a vector of percentiles.
.percentile_summary <- function(p) {
  p <- p[!is.na(p)]
  if (!length(p)) return(list(median = NA_real_, min = NA_real_, max = NA_real_))
  list(median = stats::median(p), min = min(p), max = max(p))
}

# Sorted, unique group levels, optionally restricted to scored rows.
.available_groups <- function(groups, scored = NULL) {
  g <- as.character(groups)
  if (!is.null(scored)) g <- g[!is.na(scored) & scored]
  sort(unique(stats::na.omit(g)))
}

# Two distinct existing group levels.
.check_group_levels <- function(levels, available, n = 2L) {
  if (!is.character(levels) || length(levels) != n || anyNA(levels) || anyDuplicated(levels) || !all(levels %in% available))
    stop("levels must name ", n, " distinct existing groups.", call. = FALSE)
  levels
}

.group_levels <- function(result, levels, n = 2L) {
  if (!"group" %in% names(result$scores)) stop("Score with a group column to compare groups.", call. = FALSE)
  available <- .available_groups(result$scores$group)
  if (is.null(levels)) {
    if (length(available) != n) stop("Supply levels: the data contain ", length(available), " groups.", call. = FALSE)
    return(available)
  }
  .check_group_levels(levels, available, n)
}

# Resampling units by stratum and the rows of every unit; consumes no randomness.
# `cluster` holds the values of the primary sampling unit column (one per row)
# or NULL when the rows themselves are the resampling units.
.cluster_plan <- function(n_rows, stratum, cluster) {
  unit <- if (is.null(cluster)) seq_len(n_rows) else paste(stratum, cluster, sep = "\r")
  list(units_by_stratum = lapply(split(unit, stratum), unique),
    rows_by_unit = split(seq_len(n_rows), unit))
}

# Rows of one cluster resample; must be called inside .with_seed().
.cluster_draw <- function(plan, cluster) {
  picked <- unlist(lapply(plan$units_by_stratum, function(u)
    u[sample.int(length(u), replace = TRUE)]), use.names = FALSE)
  if (is.null(cluster)) as.integer(picked) else unlist(plan$rows_by_unit[as.character(picked)], use.names = FALSE)
}
