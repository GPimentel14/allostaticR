# Specification analysis: multiverse_grid() and multiverse_analysis().

test_that("grids enumerate valid specifications and report dropped ones", {
  g <- multiverse_grid(multiverse_spec())
  expect_equal(nrow(g$choices), 6)
  expect_equal(sum(g$choices$method == "z_mean"), 2)
  expect_true(all(is.na(g$choices$risk_quantile[g$choices$method == "z_mean"])))
  expect_length(g$specs, 6); expect_equal(nrow(g$dropped), 0)
  both <- define_protocol(list(define_biomarker("a", "s", "both", "u"), define_biomarker("b", "s", "high", "u")))
  gb <- multiverse_grid(both)
  expect_equal(unique(gb$choices$method), "quantile_count")
  expect_match(gb$dropped$reason[1], "high and low")
  gc <- multiverse_grid(multiverse_spec(), method = c("quantile_count", "cutoff_count"), balance = "biomarker")
  expect_match(gc$dropped$reason, "needs a cutoff")
  expect_error(multiverse_grid(multiverse_spec(), reference = "stratified"), "strata")
  expect_error(multiverse_grid(multiverse_spec(), reference = "age"), "age column")
  expect_error(multiverse_grid(multiverse_spec(), risk_quantile = 0.4), "risk_quantile")
})

test_that("stability, agreement and specification curves are computed", {
  d <- panel_data(120, seed = 4)
  d$sbp[1] <- 190; d$hdl[1] <- 20; d$crp[1] <- 30
  d$sbp[2] <- 95; d$hdl[2] <- 90; d$crp[2] <- 0.1
  set.seed(1); d$y <- 0.05 * d$sbp + rnorm(120)
  g <- multiverse_grid(multiverse_spec())
  m <- multiverse_analysis(d, d, g, id = "id", group = "group", outcome = "y")
  expect_equal(m$summary$n_specs, 6)
  expect_equal(nrow(m$scores), 6 * 120)
  expect_equal(nrow(m$agreement), 15)
  expect_equal(m$stability$class[1:2], c("robustly high", "robustly not high"))
  expect_equal(m$stability$share_high[1:2], c(1, 0))
  expect_true(all(m$agreement$spearman > 0))
  expect_equal(nrow(m$effects), 6)
  expect_true(all(m$effects$estimate > 0))
  expect_true(all(c("stability", "agreement", "curve", "choices") %in% names(m$plots)))
  expect_silent(lapply(m$plots[c("stability", "agreement", "curve", "choices")], ggplot2::ggplot_build))
  expect_output(print(m), "6 specifications")
  d$event <- d$y > stats::median(d$y)
  logistic <- multiverse_analysis(d, d, g, outcome = "event", family = "binomial", covariates = "crp", top = 0.5)
  expect_equal(logistic$family, "binomial")
  expect_true(all(is.finite(logistic$effects$estimate)))
  expect_error(multiverse_analysis(d, d, g, outcome = "missing"), "outcome")
})

test_that("failed specifications are reported without stopping the analysis", {
  skip_if_not_installed("quantreg")
  d <- simulate_biomarkers(60, seed = 3, covariates = TRUE)
  g <- multiverse_grid(multiverse_spec(), method = "quantile_count", balance = "biomarker",
    reference = c("pooled", "stratified", "age"), strata = "sex", age = "age")
  m <- multiverse_analysis(d, d, g)
  expect_true(nrow(m$failures) >= 1)
  expect_match(m$failures$error[1], "age model")
  expect_equal(m$summary$n_specs + m$summary$n_failed, nrow(g$choices))
})
