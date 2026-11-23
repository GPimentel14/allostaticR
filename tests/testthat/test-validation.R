# Cross-validated validation of scores against outcomes: validate_scores().

test_that("cross-validated validation ranks informative scores", {
  skip_if_not_installed("survival")
  set.seed(3)
  d <- simulate_biomarkers(300, seed = 5)
  d$event <- stats::rbinom(300, 1, stats::plogis(-6 + 0.045 * d$sbp))
  d$time <- stats::rexp(300, exp(-6 + 0.03 * d$sbp)); d$age <- stats::rnorm(300, 50, 10)
  specs <- list(sbp_only = define_protocol(list(define_biomarker("sbp", "cv", "high", "mmHg"))),
    crp_only = define_protocol(list(define_biomarker("crp", "inf", "high", "mg/L"))))
  v <- validate_scores(d, specs, "event", folds = 3, covariates = "age")
  auc <- v$summary[v$summary$metric == "auc", ]
  expect_gt(auc$mean[auc$spec == "sbp_only"], auc$mean[auc$spec == "crp_only"])
  expect_true("increment_over_covariates" %in% names(v$summary))
  surv <- validate_scores(transform(d, event = 1), specs, "event", type = "survival", time = "time", folds = 3)
  expect_true(all(c("c_index", "calibration_slope") %in% surv$summary$metric))
  cont <- validate_scores(transform(d, y = sbp / 10 + stats::rnorm(300)), specs, "y", type = "continuous", folds = 3)
  expect_true("r_squared" %in% cont$summary$metric)
  expect_silent(ggplot2::ggplot_build(v$plots$performance))
})

test_that("validate_scores warns and records failing specifications", {
  set.seed(1)
  d <- data.frame(a = rnorm(30), b = rnorm(30), y = rnorm(30))
  ok <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")))
  bad <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("zzz", "s2", "high", "u")))
  expect_warning(
    v <- validate_scores(d, list(good = ok, broken = bad), outcome = "y", type = "continuous", folds = 3),
    "broken"
  )
  expect_equal(unique(v$summary$spec), "good")
  expect_true(nrow(v$failures) == 3L)
  expect_true(all(v$failures$spec == "broken"))
  expect_error(suppressWarnings(validate_scores(d, list(broken = bad), outcome = "y", type = "continuous", folds = 3)),
    "All specifications failed")
})
