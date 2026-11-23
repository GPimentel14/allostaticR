# Cohort harmonization: harmonize_cohorts().

test_that("cohort harmonization aligns distributions to the target", {
  set.seed(7)
  d <- data.frame(cohort = rep(c("main", "lab2"), c(200, 100)), a = c(rnorm(200, 50, 5), rnorm(100, 60, 10)))
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")))
  h <- harmonize_cohorts(d, spec, "cohort")
  expect_equal(h$target, "main")
  dg <- h$diagnostics[h$diagnostics$cohort == "lab2", ]
  expect_lt(dg$ks_after, dg$ks_before)
  expect_equal(dg$median_after, dg$median_target, tolerance = 0.5)
  expect_identical(h$data$a_original, d$a)
  expect_identical(h$data$a[d$cohort == "main"], d$a[d$cohort == "main"])
  ms <- harmonize_cohorts(d, spec, "cohort", method = "mean_sd")$data
  expect_equal(stats::sd(ms$a[d$cohort == "lab2"]), stats::sd(d$a[d$cohort == "main"]))
  expect_silent(ggplot2::ggplot_build(h$plots$distributions))
})

test_that("mean_sd harmonization fails clearly on zero-variance cohorts", {
  d <- data.frame(cohort = rep(c("A", "B"), each = 10), a = c(rep(1, 10), 1:10))
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")))
  expect_error(harmonize_cohorts(d, spec, "cohort", target = "B", method = "mean_sd"), "no variance")
  expect_error(harmonize_cohorts(d, spec, "cohort", target = "A", method = "mean_sd"), "no variance")
})
