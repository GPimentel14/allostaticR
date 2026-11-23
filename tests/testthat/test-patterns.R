# Risk patterns: risk_patterns().

test_that("risk patterns, co-occurrence and redundancy match hand calculations", {
  d <- data.frame(a = c(0, 10, 10, 0, 10, NA), b = c(6, 6, 5, 5, 6, 5), g = c("A", "A", "B", "B", "B", "B"))
  r <- score_load(d, fit_reference(spec = cutoff_spec("proportion", 0.5)), group = "g")
  x <- risk_patterns(r)
  expect_equal(x$n_complete, 5)
  expect_equal(x$patterns$pattern, c("a", "none", "b", "a+b"))
  expect_equal(x$patterns$n, c(2, 1, 1, 1))
  expect_equal(x$patterns$load, c(1, 0, 1, 2))
  expect_equal(x$patterns$n_systems, c(1, 0, 1, 2))
  expect_equal(x$cooccurrence$n_both, 1)
  expect_equal(x$cooccurrence$expected_both, 1.2)
  expect_equal(x$cooccurrence$phi, -1 / 6)
  expect_equal(x$redundancy$meff, 2)
  shares <- tapply(x$pattern_groups$share, x$pattern_groups$group, sum)
  expect_equal(as.vector(shares), c(1, 1))
  expect_output(print(x), "Distinct patterns: 4")
  expect_silent(lapply(x$plots, ggplot2::ggplot_build))
})

test_that("redundant indicators reduce the effective number of indicators", {
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u", 5), define_biomarker("b", "s1", "high", "u", 5),
    define_biomarker("c", "s2", "high", "u", 5)), "cutoff_count")
  d <- data.frame(a = c(1, 9, 1, 9, 9, 1), b = c(1, 9, 1, 9, 9, 1), c = c(9, 9, 1, 1, 9, 1))
  x <- risk_patterns(score_load(d, fit_reference(spec = spec)))
  expect_equal(x$redundancy$meff, 2)
  expect_true(x$cooccurrence$same_system[1]); expect_equal(x$cooccurrence$phi[1], 1)
  expect_error(risk_patterns(allostatic_load(simulate_biomarkers(20), define_protocol(list(
    define_biomarker("sbp", "s", "high", "u")), method = "z_mean"))), "binary risk flags")
})

test_that("risk_patterns disambiguates labels and plots the right members", {
  spec <- define_protocol(list(
    define_biomarker("CD4+", "imm", "high", "u", 5),
    define_biomarker("none", "met", "high", "u", 5)
  ), "cutoff_count")
  d <- data.frame("CD4+" = c(9, 1, 9, 1), "none" = c(1, 9, 9, 1), grp = c("A", "A", "B", "B"), check.names = FALSE)
  x <- risk_patterns(score_load(d, fit_reference(spec = spec), group = "grp"))
  expect_equal(sort(x$patterns$load), c(0, 1, 1, 2))
  expect_equal(sort(x$patterns$n_systems), c(0, 1, 1, 2))
  expect_equal(anyDuplicated(x$patterns$pattern), 0L)
  expect_false("key" %in% names(x$patterns))
  expect_false("rank" %in% names(x$pattern_groups))
  expect_equal(anyDuplicated(x$pattern_groups[c("group", "pattern")]), 0L)
  expect_setequal(unique(x$pattern_groups$pattern), x$patterns$pattern)
  # The member layer is layer 3 because this panel has a two-biomarker pattern.
  members <- ggplot2::ggplot_build(x$plots$matrix)$data[[3]]
  expect_equal(as.integer(table(factor(members$x, levels = seq_len(nrow(x$patterns))))), x$patterns$load)
  expect_silent(ggplot2::ggplot_build(x$plots$frequency))
})

test_that("risk_patterns keeps colliding biomarker combinations distinct", {
  spec <- define_protocol(list(
    define_biomarker("A", "s1", "high", "u", 5),
    define_biomarker("B", "s2", "high", "u", 5),
    define_biomarker("A+B", "s3", "high", "u", 5)
  ), "cutoff_count")
  d <- data.frame(A = c(1, 9, 1, 1, 9), B = c(1, 1, 9, 1, 9),
    "A+B" = c(1, 1, 1, 9, 1), check.names = FALSE)
  x <- risk_patterns(score_load(d, fit_reference(spec = spec)))
  expect_equal(nrow(x$patterns), 5L)
  expect_equal(sum(x$patterns$n), 5L)
  expect_equal(anyDuplicated(x$patterns$pattern), 0L)
  expect_silent(ggplot2::ggplot_build(x$plots$matrix))
  expect_silent(ggplot2::ggplot_build(x$plots$frequency))
})

test_that("risk_patterns rejects a group column with no values", {
  spec <- define_protocol(list(
    define_biomarker("a", "s1", "high", "u", 5),
    define_biomarker("b", "s2", "high", "u", 5)
  ), "cutoff_count")
  d <- data.frame(a = c(1, 9, 1, 9), b = c(1, 1, 9, 9), g = NA_character_)
  r <- score_load(d, fit_reference(d, spec), group = "g")
  expect_error(risk_patterns(r), "no non-missing")
})

test_that("plot_risk_map fails cleanly when there are no valid scores", {
  skip_if_not_installed("patchwork")
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")))
  ref <- fit_reference(data.frame(a = 1:10, b = 1:10), spec)
  res <- score_load(data.frame(a = c(NA_real_, NA_real_), b = c(NA_real_, NA_real_)), ref)
  expect_error(plot_risk_map(res), "No valid scores")
})
