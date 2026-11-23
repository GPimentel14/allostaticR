# Outlier detection: detect_outliers().

test_that("outliers are flagged, optionally winsorized and their impact measured", {
  skip_if_not_installed("MASS")
  set.seed(8)
  d <- data.frame(a = c(rnorm(99, 10), 40), b = rnorm(100, 5))
  d$b[50] <- 5; d$a[50] <- 13.5
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")), method = "z_mean")
  o <- detect_outliers(d, spec, winsorize = TRUE, reference = fit_reference(d, spec))
  fa <- o$flags[o$flags$biomarker == "a", ]
  expect_true(fa$outlier[100]); expect_equal(fa$direction[100], "high")
  expect_equal(o$winsorized$a[100], fa$upper_fence[100])
  expect_identical(o$winsorized$b, d$b)
  expect_true(o$impact$change[100] < 0)
  expect_true(o$multivariate$outlier[o$multivariate$row == 100])
  m <- detect_outliers(d, spec, method = "mad", multivariate = FALSE)
  expect_true(m$flags$outlier[100]); expect_null(m$winsorized)
  expect_silent(lapply(o$plots, ggplot2::ggplot_build))
  expect_output(print(o), "Outlier screen")
  expect_error(detect_outliers(d, spec, reference = fit_reference(d, spec)), "winsorize")
})
