test_that("as.data.frame methods work across allostaticR result classes", {
  dat <- simulate_biomarkers(n = 80, seed = 1)
  dat$w <- runif(nrow(dat), 0.5, 2.0)
  spec <- define_protocol(
    list(
      define_biomarker("sbp", "Cardiovascular", "high", "mmHg"),
      define_biomarker("hdl", "Metabolic", "low", "mg/dL"),
      define_biomarker("crp", "Inflammatory", "high", "mg/L")
    ),
    method = "quantile_count",
    missing = "proportion",
    min_coverage = 0.6
  )
  ref <- fit_reference(dat, spec)
  res <- score_load(dat, ref, id = "id", group = "group")

  svy <- survey_summary(res, dat, weights = "w", by = "group", B = 10L, seed = 1L)
  expect_s3_class(as.data.frame(svy), "data.frame")
  expect_equal(nrow(as.data.frame(svy)), nrow(svy$estimates))

  gap <- decompose_gap(res, levels = c("A", "B"), by = "system", B = 10L, seed = 1L)
  expect_s3_class(as.data.frame(gap), "data.frame")
  expect_s3_class(as.data.frame(gap, what = "summary"), "data.frame")

  pat <- risk_patterns(res, top = 5L)
  expect_s3_class(as.data.frame(pat), "data.frame")
  expect_s3_class(as.data.frame(pat, what = "cooccurrence"), "data.frame")
  expect_s3_class(as.data.frame(pat, what = "pattern_groups"), "data.frame")

  boot <- bootstrap_load(dat, dat, spec, B = 10L, seed = 1L)
  expect_s3_class(as.data.frame(boot), "data.frame")
  expect_s3_class(as.data.frame(boot, what = "cutoff_intervals"), "data.frame")
  expect_s3_class(as.data.frame(boot, what = "draws"), "data.frame")

  grid <- multiverse_grid(spec, method = c("quantile_count", "z_mean"), risk_quantile = 0.75)
  mv <- multiverse_analysis(dat, dat, grid, id = "id", group = "group")
  expect_s3_class(as.data.frame(mv), "data.frame")
  expect_s3_class(as.data.frame(mv, what = "stability"), "data.frame")
  expect_error(as.data.frame(mv, what = "effects"), "not available")
})
