# Bootstrap calibration: bootstrap_load().

test_that("bootstrap calibration preserves RNG, records failures and rejects cutoff protocols", {
  s <- define_protocol(list(define_biomarker("a", "s", "high", "u")))
  ref <- data.frame(a = 1:8); d <- data.frame(a = c(0, 100))
  set.seed(7); before <- .Random.seed
  b <- bootstrap_load(d, ref, s, B = 20, seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(b$draws, bootstrap_load(d, ref, s, B = 20, seed = 5)$draws)
  expect_equal(b$summary$estimate, 0.5)
  expect_equal(b$summary$lower, 0.5); expect_equal(b$summary$upper, 0.5)
  expect_equal(b$cutoff_intervals$estimate, 6.25)
  expect_silent(ggplot2::ggplot_build(b$plots$mean_interval))
  failing <- bootstrap_load(d, data.frame(a = c(1, 2, NA, NA)), s, B = 50, seed = 3)
  expect_gt(nrow(failing$failures), 0)
  expect_equal(nrow(failing$draws), 50)
  expect_equal(failing$summary$n_failed_reference, nrow(failing$failures))
  expect_error(bootstrap_load(data.frame(a = 1, b = 2), ref, cutoff_spec()), "quantile")
})

test_that("the cluster bootstrap resamples whole sampling units within strata", {
  # Two primary sampling units per stratum: resampling units mixes them,
  # resampling whole strata would keep every draw identical.
  d <- data.frame(stratum = rep(c("s1", "s2"), each = 4), psu = rep(1:2, 4),
    a = c(0, 0, 10, 20, 0, 0, 10, 20))
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u")))
  ev <- data.frame(a = c(15, 15))
  b <- bootstrap_load(ev, d, spec, B = 20, seed = 5, cluster = "psu", design_strata = "stratum")
  expect_lt(b$cutoff_intervals$lower, b$cutoff_intervals$upper)
  expect_gt(length(unique(b$draws$mean_score)), 1)
  expect_equal(b$summary$n_valid, 20)
})

test_that("the bootstrap arguments are checked before any resampling", {
  s <- define_protocol(list(define_biomarker("a", "s", "high", "u")))
  ref <- data.frame(a = 1:8); d <- data.frame(a = c(0, 100))
  expect_error(bootstrap_load(d, ref, s, B = 1), "B must be an integer >= 2.", fixed = TRUE)
  expect_error(bootstrap_load(d, ref, s, conf = 0), "conf must be in (0, 1).", fixed = TRUE)
  expect_error(bootstrap_load(d, ref, s, cluster = "missing"), "Design columns must be complete")
})
