# Sensitivity to alternative protocols: compare_scores(), sensitivity_analysis()
# and panel_influence().

test_that("sensitivity and influence have known paired normalized changes", {
  refdata <- data.frame(a = 1:4, b = 1:4)
  s <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "low", "u")))
  t <- s; t$risk_quantile <- 0.9
  r <- sensitivity_analysis(refdata, refdata, list(base = s, alternate = t))
  expect_equal(r$comparisons$n_matched, 4)
  expect_equal(r$comparisons$mean_absolute_difference, 0)
  expect_silent(ggplot2::ggplot_build(r$plots$changes))
  d <- data.frame(a = c(10, 0, NA), b = c(6, 5, 5))
  influence <- panel_influence(d, fit_reference(spec = cutoff_spec()))
  expect_equal(influence$summary$n_matched, c(2, 2))
  expect_equal(influence$summary$mean_absolute_change, c(0.5, 0.5))
  expect_equal(influence$summary$n_valid_without, c(3, 2))
  expect_equal(influence$results$a$reference$cutoffs$lower, 5)
  expect_silent(ggplot2::ggplot_build(influence$plots$influence))
  expect_error(panel_influence(d, fit_reference(spec = define_protocol(list(define_biomarker("a", "s", "high", "u", 10)), "cutoff_count"))), "two")
})

test_that("sensitivity heatmaps distinguish missing pairs from zero change", {
  d <- data.frame(a = c(0, 10, NA), b = c(6, 5, 5))
  s <- cutoff_spec(); p <- cutoff_spec("proportion", 0.5)
  x <- sensitivity_analysis(d, d, list(complete = s, partial = p))
  built <- ggplot2::ggplot_build(x$plots$changes)
  expect_equal(nrow(built$data[[1]]), 6)
  expect_true("#5b5a56" %in% built$data[[1]]$fill)
})
