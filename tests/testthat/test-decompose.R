test_that("group gaps decompose exactly into prevalence differences", {
  d <- data.frame(a = c(0, 10, 10, 10, 0), b = c(6, 6, 5, 6, 5), g = c("A", "A", "B", "B", "B"))
  r <- score_load(d, fit_reference(spec = cutoff_spec()), group = "g")
  x <- decompose_gap(r, B = 0)
  expect_equal(x$summary$mean_from, 0.5); expect_equal(x$summary$mean_to, 4 / 3)
  expect_equal(x$contributions$contribution, c(1 / 6, 2 / 3))
  expect_equal(sum(x$contributions$contribution), x$summary$difference)
  expect_equal(x$contributions$mean_from, c(0.5, 0)); expect_equal(x$contributions$mean_to, c(2 / 3, 2 / 3))
  expect_true(all(is.na(x$contributions$lower)))
  reversed <- decompose_gap(r, levels = c("B", "A"), B = 0)
  expect_equal(reversed$summary$difference, -x$summary$difference)
  expect_output(print(x), "Decomposition of B minus A")
  expect_silent(ggplot2::ggplot_build(x$plots$waterfall))
})

test_that("decomposition is exact for system-balanced and z-score scales", {
  set.seed(4)
  d <- data.frame(g = rep(c("A", "B"), each = 30), a = rnorm(60, 10), b = rnorm(60, 10), c = rnorm(60, 10))
  d$a[d$g == "B"] <- d$a[d$g == "B"] + 1
  bm <- list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s1", "high", "u"),
    define_biomarker("c", "s2", "low", "u"))
  for (spec in list(define_protocol(bm, balance = "system"), define_protocol(bm, method = "z_mean"),
                    define_protocol(bm, method = "z_mean", balance = "system"))) {
    r <- allostatic_load(d, spec, group = "g")
    x <- decompose_gap(r, B = 0)
    expect_equal(sum(x$contributions$contribution), x$summary$difference, tolerance = 1e-10)
    s <- decompose_gap(r, by = "system", B = 0)
    expect_equal(s$contributions$term, c("s1", "s2"))
    expect_equal(s$contributions$contribution[1], sum(x$contributions$contribution[1:2]), tolerance = 1e-10)
  }
})

test_that("partial panels are excluded and reported", {
  d <- data.frame(a = c(0, 10, 10, NA, 0, 10), b = c(6, 6, 5, 5, 6, 5), g = rep(c("A", "B"), each = 3))
  r <- score_load(d, fit_reference(spec = cutoff_spec("proportion", 0.5)), group = "g")
  x <- decompose_gap(r, B = 0)
  expect_equal(x$summary$n_excluded, 1)
  expect_equal(c(x$summary$n_from, x$summary$n_to), c(3, 2))
  expect_equal(sum(x$contributions$contribution), x$summary$difference)
})

test_that("bootstrap intervals are reproducible and inputs are checked", {
  set.seed(9); before <- .Random.seed
  d <- simulate_biomarkers(120)
  r <- allostatic_load(d, define_protocol(list(define_biomarker("sbp", "cv", "high", "mmHg"),
    define_biomarker("hdl", "met", "low", "mg/dL"))), group = "group")
  x <- decompose_gap(r, B = 50, seed = 3)
  expect_identical(.Random.seed, before)
  expect_identical(x$contributions, decompose_gap(r, B = 50, seed = 3)$contributions)
  expect_true(all(x$contributions$lower <= x$contributions$contribution))
  expect_true(x$summary$lower <= x$summary$difference && x$summary$difference <= x$summary$upper)
  expect_error(decompose_gap(allostatic_load(d, r$reference$spec)), "group column")
  three <- d; three$group <- rep(c("A", "B", "C"), 40)
  r3 <- allostatic_load(three, r$reference$spec, group = "group")
  expect_error(decompose_gap(r3), "Supply levels")
  expect_equal(decompose_gap(r3, levels = c("A", "C"), B = 0)$summary$to, "C")
  expect_error(decompose_gap(r3, levels = c("A", "Z")), "levels")
})
