# Internal helpers: regression tests for the shared implementation details.
# The helpers are not exported; they are reached with allostaticR:::.

test_that(".check_integer preserves custom messages and enforces bounds", {
  expect_error(allostaticR:::.check_integer(1L, "B", 2L, text = "B must be an integer >= 2."),
    "B must be an integer >= 2.", fixed = TRUE)
  expect_error(allostaticR:::.check_integer(0L, "n", 1L, text = "n must be a positive integer."),
    "n must be a positive integer.", fixed = TRUE)
  expect_error(allostaticR:::.check_integer(1.5, "B", 2L), "B must be an integer from 2 upwards.", fixed = TRUE)
  expect_error(allostaticR:::.check_integer(5L, "n", 1L, 3L), "n must be an integer from 1 to 3.", fixed = TRUE)
  expect_error(allostaticR:::.check_integer(NA_integer_, "n", 1L), "must be an integer")
  expect_error(allostaticR:::.check_integer(c(1L, 2L), "n", 1L), "must be an integer")
  expect_error(allostaticR:::.check_integer("2", "B", 2L), "must be an integer")
  expect_error(allostaticR:::.check_integer(Inf, "B", 2L), "must be an integer")
  expect_no_error(allostaticR:::.check_integer(2L, "B", 2L))
  # the public entry points keep their original wording
  expect_error(simulate_biomarkers(0), "n must be a positive integer.", fixed = TRUE)
  s <- define_protocol(list(define_biomarker("a", "s", "high", "u")))
  expect_error(bootstrap_load(data.frame(a = 1:8), data.frame(a = 1:8), s, B = 1),
    "B must be an integer >= 2.", fixed = TRUE)
})

test_that(".check_seed accepts zero and rejects invalid seeds", {
  expect_error(allostaticR:::.check_seed(-1), "seed must be a non-negative integer.", fixed = TRUE)
  expect_error(allostaticR:::.check_seed(1.5), "seed must be a non-negative integer.", fixed = TRUE)
  expect_error(allostaticR:::.check_seed(NA), "seed must be a non-negative integer.", fixed = TRUE)
  expect_error(allostaticR:::.check_seed("1"), "seed must be a non-negative integer.", fixed = TRUE)
  expect_error(allostaticR:::.check_seed(c(1, 2)), "seed must be a non-negative integer.", fixed = TRUE)
  expect_error(allostaticR:::.check_seed(Inf), "seed must be a non-negative integer.", fixed = TRUE)
  expect_no_error(allostaticR:::.check_seed(0))
  expect_no_error(allostaticR:::.check_seed(1L))
})

test_that(".check_data_frame rejects empty data and duplicated names", {
  expect_error(allostaticR:::.check_data_frame(list(a = 1)), "data must be a non-empty data frame.", fixed = TRUE)
  expect_error(allostaticR:::.check_data_frame(data.frame()), "data must be a non-empty data frame.", fixed = TRUE)
  expect_error(allostaticR:::.check_data_frame(data.frame(), "frame"),
    "frame must be a non-empty data frame.", fixed = TRUE)
  duplicated <- data.frame(a = 1, b = 2); names(duplicated) <- c("a", "a")
  expect_error(allostaticR:::.check_data_frame(duplicated), "Column names must be unique.", fixed = TRUE)
  expect_no_error(allostaticR:::.check_data_frame(data.frame(a = 1)))
})

test_that(".invalid_medication_columns reports absent, non-logical and incomplete columns", {
  d <- data.frame(a = 1, treated = c(TRUE, FALSE))
  expect_identical(allostaticR:::.invalid_medication_columns(d, character()), character())
  expect_identical(allostaticR:::.invalid_medication_columns(d, c(a = "treated")), character())
  expect_identical(allostaticR:::.invalid_medication_columns(d, c(a = "missing")), "missing")
  d$flag <- c(TRUE, NA)
  expect_identical(allostaticR:::.invalid_medication_columns(d, c(a = "flag")), "flag")
  d$num <- c(1, 0)
  expect_identical(allostaticR:::.invalid_medication_columns(d, c(a = "num")), "num")
  expect_identical(allostaticR:::.invalid_medication_columns(d, c(a = "treated", "missing")), "missing")
})

test_that(".risk_flag applies one rule for every direction and inclusivity", {
  p <- list(lower = 5, upper = 10)
  x <- c(4, 5, 7, 10, 11)
  expect_identical(allostaticR:::.risk_flag(x, p, "low", TRUE), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(allostaticR:::.risk_flag(x, p, "low", FALSE), c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(allostaticR:::.risk_flag(x, p, "high", TRUE), c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(allostaticR:::.risk_flag(x, p, "high", FALSE), c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(allostaticR:::.risk_flag(x, p, "both", TRUE), c(TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_identical(allostaticR:::.risk_flag(x, p, "both", FALSE), c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(allostaticR:::.risk_flag(c(1, NA), p, "high", TRUE), c(FALSE, NA))
  # the flags of score_load() are the flags of the helper
  spec <- cutoff_spec()
  ref <- fit_reference(spec = spec)
  cut <- reference_cutoffs(ref)
  d <- data.frame(a = c(11, 9, 10, NA), b = c(4, 6, 5, 5))
  components <- as.data.frame(score_load(d, ref), what = "components")
  a <- components[components$biomarker == "a", ]
  b <- components[components$biomarker == "b", ]
  expect_identical(a$risk, as.integer(allostaticR:::.risk_flag(d$a,
    list(lower = cut$lower[1], upper = cut$upper[1]), "high", spec$inclusive)))
  expect_identical(b$risk, as.integer(allostaticR:::.risk_flag(d$b,
    list(lower = cut$lower[2], upper = cut$upper[2]), "low", spec$inclusive)))
})

test_that(".cluster_plan and .cluster_draw reproduce the previous inline resampling", {
  d <- data.frame(stratum = rep(c("s1", "s2"), each = 6), psu = rep(1:3, 4))
  n <- nrow(d)
  stratum <- as.character(d$stratum)
  plan <- allostaticR:::.cluster_plan(n, stratum, d$psu)
  expect_named(plan, c("units_by_stratum", "rows_by_unit"))
  expect_length(plan$units_by_stratum, 2L)
  expect_length(plan$rows_by_unit, 6L)
  legacy <- function(seed) {
    allostaticR:::.with_seed(seed, function() {
      unit <- paste(stratum, d$psu, sep = "\r")
      units_by_stratum <- lapply(split(unit, stratum), unique)
      rows_by_unit <- split(seq_len(n), unit)
      picked <- unlist(lapply(units_by_stratum, function(u) u[sample.int(length(u), replace = TRUE)]),
        use.names = FALSE)
      unlist(rows_by_unit[as.character(picked)], use.names = FALSE)
    })
  }
  set.seed(101); before <- .Random.seed
  rows <- allostaticR:::.with_seed(42, function() allostaticR:::.cluster_draw(plan, "psu"))
  expect_identical(.Random.seed, before)
  expect_identical(rows, legacy(42))
  expect_equal(rows, c(1, 4, 1, 4, 1, 4, 7, 10, 8, 11, 8, 11))
  # without clusters the rows themselves are the resampling units
  row_plan <- allostaticR:::.cluster_plan(n, stratum, NULL)
  row_rows <- allostaticR:::.with_seed(7, function() allostaticR:::.cluster_draw(row_plan, NULL))
  expect_identical(row_rows, allostaticR:::.with_seed(7, function() {
    unit <- seq_len(n)
    units_by_stratum <- lapply(split(unit, stratum), unique)
    picked <- unlist(lapply(units_by_stratum, function(u) u[sample.int(length(u), replace = TRUE)]),
      use.names = FALSE)
    as.integer(picked)
  }))
})

test_that(".midrank_percentile shares mid-ranks between ties", {
  expect_equal(allostaticR:::.midrank_percentile(c(1, 2, 2, 4)), c(0.125, 0.5, 0.5, 0.875))
  expect_equal(allostaticR:::.midrank_percentile(c(1, 2, 3)), (rank(c(1, 2, 3)) - 0.5) / 3)
  expect_equal(allostaticR:::.midrank_percentile(c(5, 5, 5)), rep(0.5, 3))
  expect_true(is.na(allostaticR:::.midrank_percentile(numeric())))
  expect_false(is.nan(allostaticR:::.midrank_percentile(numeric())))
})

test_that(".paired_scores matches the comparison of two score vectors", {
  x <- c(1, 2, 3, NA); y <- c(2, 2, 4, 5)
  p <- allostaticR:::.paired_scores(x, y)
  expect_equal(p$n_matched, 3)
  expect_equal(p$mean_difference, mean(y[1:3] - x[1:3]))
  expect_equal(p$mean_absolute_difference, mean(abs(y[1:3] - x[1:3])))
  expect_equal(p$spearman, stats::cor(x[1:3], y[1:3], method = "spearman"))
  expect_true(is.na(allostaticR:::.paired_scores(c(1, 1), c(1, 1))$spearman))
  # the same numbers as compare_scores()
  ref <- fit_reference(spec = cutoff_spec())
  d <- data.frame(id = c("a", "b", "c"), a = c(10, 0, 11), b = c(5, 6, 4))
  r1 <- score_load(d, ref, id = "id")
  r2 <- score_load(d, ref, id = "id")
  cmp <- compare_scores(list(first = r1, second = r2))
  paired <- allostaticR:::.paired_scores(r1$scores$score, r2$scores$score)
  expect_equal(cmp$n_matched, paired$n_matched)
  expect_equal(cmp$mean_difference_b_minus_a, paired$mean_difference)
  expect_equal(cmp$mean_absolute_difference, paired$mean_absolute_difference)
  expect_equal(cmp$spearman, paired$spearman)
})

test_that(".palette_for recycles colours for any length", {
  cols <- c("#111111", "#222222")
  expect_identical(allostaticR:::.palette_for(3, cols), c("#111111", "#222222", "#111111"))
  expect_identical(allostaticR:::.palette_for(2, cols), cols)
  expect_identical(allostaticR:::.palette_for(0, cols), character())
  expect_length(allostaticR:::.palette_for(5, cols), 5L)
  expect_identical(allostaticR:::.palette_for(2), allostaticR:::.al_palette$categorical[1:2])
})

test_that(".random_subset is reproducible and restores the caller's random state", {
  set.seed(11); before <- .Random.seed
  s1 <- allostaticR:::.random_subset(1:10, 3, seed = 7)
  expect_identical(.Random.seed, before)
  expect_identical(s1, allostaticR:::.random_subset(1:10, 3, seed = 7))
  expect_identical(s1, c(3L, 7L, 10L))
  expect_identical(s1, sort(s1))
  expect_length(s1, 3L)
  # no randomness is consumed when the whole vector is taken
  set.seed(11); before <- .Random.seed
  expect_identical(allostaticR:::.random_subset(1:10, 10, seed = 7), 1:10)
  expect_identical(allostaticR:::.random_subset(1:10, 20, seed = 7), 1:10)
  expect_identical(.Random.seed, before)
})

test_that(".print_table prints without row names and with fixed digits", {
  d <- data.frame(a = 1.23456, b = 2)
  expect_identical(capture.output(allostaticR:::.print_table(d)),
    capture.output(print(d, row.names = FALSE, digits = 3)))
  expect_identical(capture.output(allostaticR:::.print_table(d, digits = 1)),
    capture.output(print(d, row.names = FALSE, digits = 1)))
  expect_false(any(grepl("^1 ", capture.output(allostaticR:::.print_table(d)))))
})

test_that(".percentile_summary reports median, minimum and maximum", {
  s <- allostaticR:::.percentile_summary(c(0.1, 0.5, 0.9))
  expect_equal(s$median, 0.5); expect_equal(s$min, 0.1); expect_equal(s$max, 0.9)
  expect_equal(allostaticR:::.percentile_summary(c(NA, 0.2, 0.4))$median, 0.3)
  for (empty in list(c(NA_real_, NA_real_), numeric())) {
    summary <- allostaticR:::.percentile_summary(empty)
    expect_true(is.na(summary$median)); expect_true(is.na(summary$min)); expect_true(is.na(summary$max))
  }
})

test_that(".available_groups and .check_group_levels share one group policy", {
  expect_identical(allostaticR:::.available_groups(c("b", "a", "b", NA)), c("a", "b"))
  expect_identical(allostaticR:::.available_groups(c("b", "a", NA), scored = c(TRUE, FALSE, TRUE)), "b")
  expect_identical(allostaticR:::.check_group_levels(c("a", "b"), c("a", "b")), c("a", "b"))
  expect_error(allostaticR:::.check_group_levels(c("a", "z"), c("a", "b")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
  expect_error(allostaticR:::.check_group_levels(c("a", "a"), c("a", "b")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
  expect_error(allostaticR:::.check_group_levels("a", c("a", "b")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
  expect_error(allostaticR:::.check_group_levels(c("a", NA), c("a", "b")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
})

test_that(".group_levels keeps its strict two-group policy", {
  d <- data.frame(a = c(0, 10, 0, 10), b = c(5, 6, 4, 5), g = c("A", "A", "B", "B"))
  r <- score_load(d, fit_reference(spec = cutoff_spec()), group = "g")
  expect_identical(allostaticR:::.group_levels(r, NULL), c("A", "B"))
  expect_identical(allostaticR:::.group_levels(r, c("B", "A")), c("B", "A"))
  expect_error(allostaticR:::.group_levels(r, NULL, n = 3L),
    "Supply levels: the data contain 2 groups.", fixed = TRUE)
  expect_error(allostaticR:::.group_levels(r, c("A", "Z")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
})

test_that(".biomarker_field keeps the protocol order", {
  spec <- panel_spec()
  expect_identical(allostaticR:::.biomarker_field(spec, "name"), c("sbp", "hdl", "crp"))
  expect_identical(allostaticR:::.biomarker_field(spec, "system"),
    c("cardiovascular", "metabolic", "inflammatory"))
  expect_identical(allostaticR:::.biomarker_field(spec, "label"), c("sbp", "hdl", "crp"))
  expect_identical(allostaticR:::.biomarker_field(spec, "transform"), rep("identity", 3))
  expect_identical(allostaticR:::.biomarker_field(spec, "direction"),
    c("high", "low", "high"))
})

test_that(".per_biomarker stacks one row per biomarker without adding columns", {
  spec <- panel_spec()
  out <- allostaticR:::.per_biomarker(spec, function(b) data.frame(biomarker = b$name, value = 1))
  expect_identical(out$biomarker, c("sbp", "hdl", "crp"))
  expect_named(out, c("biomarker", "value"))
  expect_identical(rownames(out), as.character(1:3))
  expect_false(identical(rownames(out), c("sbp", "hdl", "crp")))
  expect_equal(nrow(out), 3L)
  # the public summaries keep their own columns
  q <- check_data_quality(data.frame(sbp = 1, hdl = 2, crp = 3), spec)
  expect_named(q$summary, c("biomarker", "present", "numeric", "n", "n_missing", "n_nonfinite",
    "n_outside_limits", "expected_unit", "declared_unit", "unit_status", "missing_fraction"))
})

test_that(".finish_plot sets labels without erasing an existing title", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::labs(title = "Original title")
  out <- allostaticR:::.finish_plot(p, title = "New title", subtitle = "Sub", caption = "A caption")
  expect_equal(out$labels$title, "New title")
  expect_equal(out$labels$subtitle, "Sub")
  expect_match(out$labels$caption, "A caption")
  kept <- allostaticR:::.finish_plot(p)
  expect_equal(kept$labels$title, "Original title")
  expect_null(kept$labels$subtitle)
  expect_null(kept$labels$caption)
  # the package theme is applied
  expect_equal(kept$theme$legend.position, "bottom")
  expect_equal(kept$theme$plot.title$face, "bold")
})

test_that(".heatmap_text_colour switches to white above the threshold", {
  primary <- allostaticR:::.al_palette$primary
  expect_identical(allostaticR:::.heatmap_text_colour(c(0.9, 0.1, NA), 0.55),
    c("white", primary, primary))
  expect_identical(allostaticR:::.heatmap_text_colour(c(0.9), 0.75), "white")
  expect_identical(allostaticR:::.heatmap_text_colour(c(0.7), 0.75), primary)
  expect_identical(allostaticR:::.heatmap_text_colour(c(0.75), 0.75), primary)
})

test_that(".heatmap_layer builds tiles and centred labels", {
  d <- data.frame(col = c("x", "y"), row = c("r1", "r2"), value = c(0.9, 0.1), label = c("0.9", "0.1"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = col, y = row)) +
    allostaticR:::.heatmap_layer(d, "col", "row", "value", label = "label", threshold = 0.55)
  built <- ggplot2::ggplot_build(p)
  expect_length(built$data, 2L)
  expect_equal(nrow(built$data[[1]]), 2L)
  expect_true("fill" %in% names(built$data[[1]]))
  expect_true("colour" %in% names(built$data[[2]]))
  expect_false("fill" %in% names(built$data[[2]]))
  expect_equal(built$data[[2]]$colour,
    c("white", allostaticR:::.al_palette$primary))
  expect_equal(built$data[[2]]$label, c("0.9", "0.1"))
  # without a threshold every label uses the primary colour
  plain <- ggplot2::ggplot(d, ggplot2::aes(x = col, y = row)) +
    allostaticR:::.heatmap_layer(d, "col", "row", "value", label = "label")
  expect_equal(ggplot2::ggplot_build(plain)$data[[2]]$colour,
    rep(allostaticR:::.al_palette$primary, 2))
})

test_that(".percent matches the previous inline percentage labels", {
  x <- c(0, 0.005, 0.125, 0.5, 1)
  expect_identical(allostaticR:::.percent(x), paste0(round(x * 100), "%"))
  expect_identical(allostaticR:::.percent(x), c("0%", "0%", "12%", "50%", "100%"))
  expect_identical(allostaticR:::.percent(c(0.2, 0.25)), c("20%", "25%"))
})

test_that(".mean_or_na never lets infinities reach the mean", {
  expect_equal(allostaticR:::.mean_or_na(c(1, 2, 3)), 2)
  expect_equal(allostaticR:::.mean_or_na(c(1, 2, NA)), 1.5)
  expect_equal(allostaticR:::.mean_or_na(c(1, 2, Inf)), 1.5)
  expect_true(is.na(allostaticR:::.mean_or_na(c(Inf, -Inf))))
  expect_true(is.na(allostaticR:::.mean_or_na(c(NA, NA))))
  expect_true(is.na(allostaticR:::.mean_or_na(numeric())))
})

test_that("the overwrite argument reports one message in every exporter", {
  expect_error(save_protocol(fit_reference(spec = cutoff_spec()), tempfile(), overwrite = "yes"),
    "overwrite must be TRUE or FALSE.", fixed = TRUE)
  expect_error(allostaticR:::.check_flag(NA, "overwrite"),
    "overwrite must be TRUE or FALSE.", fixed = TRUE)
})

test_that(".scale_label keeps its own label family apart from the longitudinal ones", {
  expect_identical(allostaticR:::.scale_label("count"), "Allostatic load (count)")
  expect_identical(allostaticR:::.scale_label("proportion"), "At-risk proportion")
  expect_identical(allostaticR:::.scale_label("system_proportion"), "System-balanced risk proportion")
  expect_identical(allostaticR:::.scale_label("z"), "Mean risk-oriented z-score")
  expect_identical(allostaticR:::.scale_label("system_z"), "System-balanced mean z-score")
  expect_identical(allostaticR:::.scale_label("something else"), "something else")
  # the longitudinal graphics use another family and must not be unified with it
  d <- data.frame(person = c("A", "A", "B", "B"), time = c(0, 1, 0, 1),
    a = c(0, 10, 0, 10), b = c(5, 6, 5, 6))
  tr <- score_trajectory(d, fit_reference(spec = cutoff_spec()), "person", "time", subjects = "A")
  expect_equal(tr$plots$cohort$labels$y, "Mean score (count)")
  expect_equal(tr$plots$individual$labels$y, "Score (count)")
  expect_false(identical(tr$plots$cohort$labels$y, allostaticR:::.scale_label("count")))
  expect_false(identical(tr$plots$individual$labels$y, allostaticR:::.scale_label("count")))
})

test_that("captions are wrapped in the quality and outlier graphics", {
  q <- check_data_quality(data.frame(a = 1, b = 2), cutoff_spec())
  expect_identical(q$plots$missingness$labels$caption, allostaticR:::.wrap_caption(
    "Source: supplied data. Absent columns counted as fully unavailable; no imputation."))
  o <- detect_outliers(data.frame(a = stats::rnorm(20), b = stats::rnorm(20)),
    define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u"))))
  expect_match(o$plots$univariate$labels$caption, "Values are flagged, not removed")
  expect_identical(o$plots$univariate$labels$caption, allostaticR:::.wrap_caption(
    "Tukey fence rule with k = 3 on the protocol transform scale. Values are flagged, not removed.", 100))
})

test_that("the palette keeps the exact colours of the published figures", {
  p <- allostaticR:::.al_palette
  expect_identical(p$groups[["ink"]], "#262626")
  expect_identical(p$groups[["mid"]], "#7f7f7f")
  expect_identical(p$groups[["light"]], "#bfbfbf")
  expect_identical(p$light_grey, "#a6a6a6")
  expect_identical(p$light_blue, "#d9e2ec")
  expect_identical(p$mid_grey, "#595959")
  expect_identical(p$grey, "#8c8c8c")
  expect_identical(p$light_red, "#efd9d9")
  expect_identical(p$primary, "#0b0b0b")
  expect_identical(p$unavailable, "#5b5a56")
})
