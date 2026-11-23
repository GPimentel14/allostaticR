# Longitudinal analysis: score_trajectory(), visit_changes() and model_trajectory().

test_that("mixed models recover trajectory slopes", {
  skip_if_not_installed("nlme")
  set.seed(4)
  ids <- rep(1:60, each = 4); time <- rep(0:3, 60); g <- rep(c("A", "B"), each = 120)
  person <- rep(stats::rnorm(60, 0, 3), each = 4)
  d <- data.frame(id = ids, time = time, g = g, a = 120 + person + ifelse(g == "B", 4, 1) * time + stats::rnorm(240, 0, 1))
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")), method = "z_mean")
  ref <- fit_reference(d[d$time == 0, ], spec)
  m <- model_trajectory(score_trajectory(d, ref, "id", "time", group = "g"))
  sd0 <- reference_cutoffs(ref)$scale
  expect_equal(m$fixed$estimate[m$fixed$term == "time"], 1 / sd0, tolerance = 0.15)
  expect_equal(m$fixed$estimate[m$fixed$term == "time:groupB"], 3 / sd0, tolerance = 0.15)
  expect_silent(ggplot2::ggplot_build(m$plots$trajectory))
  expect_output(print(m), "Mixed-effects")
})

test_that("longitudinal changes use actual endpoints and reject duplicate visits", {
  d <- data.frame(person = c("A", "A", "B", "B", "C"),
    time = as.Date("2026-01-01") + c(0, 30, 0, 30, 0),
    a = c(0, 10, NA, 10, 0), b = c(6, 5, 6, 5, 6))
  tr <- score_trajectory(d, fit_reference(spec = cutoff_spec()), "person", "time", subjects = "A")
  changes <- visit_changes(tr)
  expect_equal(changes$change, c(2, NA, NA))
  expect_equal(changes$status, c("ok", "missing_endpoint", "single_visit"))
  expect_equal(changes$elapsed, c(30, 30, 0))
  expect_equal(tr$summary$n_valid, c(2, 2))
  expect_silent(lapply(tr$plots, ggplot2::ggplot_build))
  expect_error(score_trajectory(rbind(d, d[1, ]), tr$reference, "person", "time"), "Duplicate")
  d$time[1] <- NA
  expect_error(score_trajectory(d, tr$reference, "person", "time"), "finite")
})

test_that("score_trajectory rejects a group column with no values", {
  d <- data.frame(person = c("A", "A", "B", "B"), time = c(0, 1, 0, 1),
    a = c(0, 10, 10, 0), b = c(5, 6, 4, 4), g = NA_character_)
  ref <- fit_reference(d[d$time == 0, ], cutoff_spec())
  expect_error(score_trajectory(d, ref, "person", "time", group = "g"), "no non-missing")
})

test_that("model_trajectory warns and updates label when falling back to random intercept", {
  skip_if_not_installed("nlme")
  set.seed(4)
  d <- data.frame(id = rep(1:10, each = 3), time = rep(0:2, 10), a = rnorm(30, 10, 2))
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")), method = "z_mean")
  tr <- score_trajectory(d, fit_reference(d[d$time == 0, ], spec), "id", "time")
  real_lme <- nlme::lme
  testthat::with_mocked_bindings(
    lme = function(fixed, data, random, ...) {
      if (identical(deparse(random), "~t | id")) stop("singular precision matrix")
      real_lme(fixed = fixed, data = data, random = random, ...)
    },
    .package = "nlme",
    {
      expect_warning(m <- model_trajectory(tr, random = "slope"), "random intercepts only")
      expect_equal(m$random_structure, "intercept")
      expect_match(m$plots$trajectory$labels$subtitle, "random intercepts;")
    }
  )
})
