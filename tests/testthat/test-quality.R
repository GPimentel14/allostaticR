# Data quality: check_data_quality() and convert_units().

test_that("quality and conversions expose rather than repair issues", {
  d <- data.frame(id = c(1, 1, 3), a = c(0, 10, NA), b = c(5, Inf, 6))
  q <- check_data_quality(d, cutoff_spec(), limits = list(a = c(1, 9)), units = c(a = "wrong", b = "u"), id = "id")
  expect_equal(q$summary$n_outside_limits[1], 2)
  expect_true(all(c("unit_mismatch", "nonfinite_values", "duplicate_identifiers") %in% q$issues))
  expect_silent(ggplot2::ggplot_build(q$plots$missingness))
  converted <- convert_units(data.frame(a = c(1, NA)), "a", "old", "u", factor = 10)
  expect_equal(converted$a, c(10, NA)); expect_equal(attr(converted, "al_units"), c(a = "u"))
  expect_error(convert_units(converted, "a", "other", "u", 1), "mismatch")
  wrong <- convert_units(data.frame(a = 1, b = 6), "a", "u", "other", 2)
  expect_error(score_load(wrong, fit_reference(spec = cutoff_spec())), "unit mismatch")
  r <- allostatic_load(data.frame(a = c(0, 10), b = c(5, 6)), cutoff_spec())
  expect_silent(ggplot2::ggplot_build(plot_load(r, "cutoffs")))
})

test_that("the quality summary keeps its public shape and never repairs data", {
  d <- data.frame(a = c(0, 10, NA), b = c(5, Inf, 6))
  q <- check_data_quality(d, cutoff_spec())
  expect_named(q$summary, c("biomarker", "present", "numeric", "n", "n_missing", "n_nonfinite",
    "n_outside_limits", "expected_unit", "declared_unit", "unit_status", "missing_fraction"))
  expect_equal(q$summary$biomarker, c("a", "b"))
  expect_equal(q$summary$n, c(3, 3))
  expect_equal(q$summary$n_missing, c(1, 0))
  expect_equal(q$summary$n_nonfinite, c(0, 1))
  expect_true("nonfinite_values" %in% q$issues)
  # the input data frame is not modified
  expect_identical(d, data.frame(a = c(0, 10, NA), b = c(5, Inf, 6)))
})
