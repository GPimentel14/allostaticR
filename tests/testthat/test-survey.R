# Survey analysis: survey_summary().

test_that("cluster bootstrap and survey summaries follow the design", {
  set.seed(2)
  d <- data.frame(psu = rep(1:20, each = 5), stratum = rep(1:2, each = 50), w = runif(100, 0.5, 2),
    a = rnorm(100), b = rnorm(100), g = rep(c("x", "y"), 50))
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")), weights = "w")
  boot <- bootstrap_load(d, d, spec, B = 10, cluster = "psu", design_strata = "stratum")
  expect_equal(boot$summary$n_valid, 10)
  r <- score_load(d, fit_reference(d, spec), group = "g")
  s <- survey_summary(r, d, "w", design_strata = "stratum", cluster = "psu", by = "g", B = 20)
  m <- s$estimates[s$estimates$target == "mean score" & s$estimates$scope == "total", ]
  expect_equal(m$estimate, stats::weighted.mean(r$scores$score, d$w))
  expect_true(m$lower <= m$estimate && m$estimate <= m$upper)
  expect_equal(nrow(s$estimates), 3 * 3)
  expect_silent(ggplot2::ggplot_build(s$plots$prevalence))
  expect_output(print(s), "Survey-weighted")
})

test_that("the survey bootstrap resamples primary sampling units within strata", {
  set.seed(2)
  d <- data.frame(psu = rep(1:2, 4), stratum = rep(c("s1", "s2"), each = 4), w = c(1, 2, 3, 4, 1, 2, 3, 4),
    a = c(1, 2, 3, 4, 5, 6, 7, 8), b = c(2, 3, 4, 5, 6, 7, 8, 9), g = rep(c("x", "y"), 4))
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")))
  r <- score_load(d, fit_reference(d, spec), group = "g")
  s <- survey_summary(r, d, "w", design_strata = "stratum", cluster = "psu", by = "g", B = 20, seed = 1L)
  # the same draws as the inline resampling of primary sampling units
  score <- r$scores$score; w <- d$w; stratum <- as.character(d$stratum)
  unit <- paste(stratum, d$psu, sep = "\r")
  units_by_stratum <- lapply(split(unit, stratum), unique)
  rows_by_unit <- split(seq_len(nrow(d)), unit)
  legacy <- allostaticR:::.with_seed(1L, function() vapply(seq_len(20), function(i) {
    picked <- unlist(lapply(units_by_stratum, function(u) u[sample.int(length(u), replace = TRUE)]),
      use.names = FALSE)
    rows <- unlist(rows_by_unit[as.character(picked)], use.names = FALSE)
    ok <- rows[!is.na(score[rows])]
    stats::weighted.mean(score[ok], w[ok])
  }, numeric(1)))
  expected <- allostaticR:::.percentile_interval(legacy, 0.95)
  m <- s$estimates[s$estimates$target == "mean score" & s$estimates$scope == "total", ]
  expect_equal(c(m$lower, m$upper), expected)
})

test_that("the survey design columns are checked before any estimate", {
  d <- data.frame(a = 1:4, b = 2:5, w = c(1, 1, 1, 1), g = rep(c("x", "y"), 2))
  r <- score_load(d, fit_reference(d, define_protocol(list(define_biomarker("a", "s1", "high", "u"),
    define_biomarker("b", "s2", "high", "u")))), group = "g")
  expect_error(survey_summary(r, d, "missing"), "Missing column")
  expect_error(survey_summary(r, d, "w", design_strata = "missing"), "Missing column")
  expect_error(survey_summary(r, d, "w", cluster = "missing"), "Missing column")
  bad <- d; bad$w[1] <- 0
  expect_error(survey_summary(r, bad, "w"), "Weights must be positive")
  expect_error(survey_summary(r, d, "w", B = -1), "B must be an integer")
})

test_that("survey_summary keeps a group named 'All participants' distinct from total", {
  d <- data.frame(a = c(1, 2, 10, 10), b = c(1, 2, 10, 10), w = c(1, 1, 1, 1),
    g = c("All participants", "All participants", "Other", "Other"))
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u", 5), define_biomarker("b", "s2", "high", "u", 5)), "cutoff_count")
  r <- score_load(d, fit_reference(spec = spec), group = "g")
  s <- survey_summary(r, d, "w", by = "g", B = 0)
  ms <- s$estimates[s$estimates$target == "mean score", ]
  expect_equal(nrow(ms), 3L)
  expect_equal(ms$scope, c("total", "group", "group"))
  expect_equal(ms$estimate, c(1, 0, 2))
  expect_equal(ms$group, c("All participants", "All participants", "Other"))
  expect_equal(unique(s$plots$prevalence$data$series),
    c("Total", "Group: All participants", "Group: Other"))
  expect_silent(ggplot2::ggplot_build(s$plots$prevalence))
})
