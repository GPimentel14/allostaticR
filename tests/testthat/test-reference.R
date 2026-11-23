# Reference fitting and comparison: fit_reference(), reference_cutoffs(),
# compare_references() and the protocol passports (save_protocol/read_protocol).

test_that("survey weights shift reference quantiles and means", {
  b <- list(define_biomarker("a", "s", "high", "u"))
  d <- data.frame(a = 1:4, w = c(1, 1, 1, 5))
  equal <- fit_reference(transform(d, w = 1), define_protocol(b, weights = "w"))
  expect_equal(reference_cutoffs(equal)$upper, unname(stats::quantile(1:4, 0.75, type = 5)))
  heavy <- fit_reference(d, define_protocol(b, weights = "w"))
  expect_gt(reference_cutoffs(heavy)$upper, reference_cutoffs(equal)$upper)
  z <- fit_reference(d, define_protocol(b, method = "z_mean", weights = "w"))
  expect_equal(reference_cutoffs(z)$center, stats::weighted.mean(1:4, d$w))
  expect_error(fit_reference(d[1], define_protocol(b, weights = "w")), "weights column")
  expect_error(fit_reference(transform(d, w = 0), define_protocol(b, weights = "w")), "positive")
  expect_equal(score_load(data.frame(a = 4), heavy)$scores$score, 1)
})

test_that("protocol passports preserve references and reject unsupported schema", {
  ref <- fit_reference(spec = cutoff_spec())
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  save_protocol(ref, f)
  expect_error(save_protocol(ref, f), "exists")
  saveRDS(list(schema_version = 99L, protocol = ref), f)
  expect_error(read_protocol(f), "schema")
  bad <- ref; bad$cutoffs$upper[1] <- Inf
  expect_error(save_protocol(bad, f, overwrite = TRUE), "threshold")
  bad <- cutoff_spec(); bad$biomarkers[[1]]$cutoff <- "bad"
  expect_error(allostaticR:::validate_protocol(bad), "cutoff")
})

test_that("comparability detects panel, units, medications and thresholds", {
  ref <- fit_reference(spec = cutoff_spec()); s <- cutoff_spec(); s$biomarkers[[1]]$cutoff <- 11
  changed <- fit_reference(spec = s)
  x <- compare_references(list(original = ref, identical = ref, changed = changed))
  expect_true(x$compatible[1]); expect_false(x$compatible[2])
  expect_true(x$algorithm_equal[2]); expect_false(x$thresholds_equal[2])
  d <- data.frame(id = c("a", "b", "c"), a = c(10, 0, 11), b = c(5, 6, 4))
  r <- score_load(d, ref, id = "id")
  reversed <- score_load(d[3:1, ], ref, id = "id")
  cmp <- compare_scores(list(first = r, reordered = reversed))
  expect_equal(cmp$spearman, 1); expect_equal(cmp$mean_absolute_difference, 0)
  prop <- score_load(d, fit_reference(spec = cutoff_spec("proportion", 1)), id = "id")
  expect_true(is.na(compare_scores(list(count = r, proportion = prop))$mean_difference_b_minus_a))
  expect_error(compare_scores(list(first = r, subset = score_load(d[1:2, ], ref, id = "id"))), "identifier")
})

test_that("stratified quantile references are fitted and applied per stratum", {
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")), strata = "sex")
  ref <- fit_reference(data.frame(sex = rep(c("F", "M"), each = 4), a = c(1:4, 11:14)), spec)
  cut <- reference_cutoffs(ref)
  expect_equal(cut$stratum, c("F", "M")); expect_equal(cut$upper, c(3.25, 13.25))
  r <- score_load(data.frame(sex = c("F", "M", "F", "M"), a = c(4, 12, 3, 14)), ref)
  expect_equal(r$scores$score, c(1, 0, 0, 1))
  expect_equal(r$scores$stratum, c("F", "M", "F", "M"))
  expect_error(score_load(data.frame(sex = "X", a = 1), ref), "absent from the reference")
  expect_error(score_load(data.frame(sex = NA, a = 1), ref), "complete")
  expect_error(score_load(data.frame(a = 1), ref), "strata column")
  expect_error(fit_reference(data.frame(sex = c("F", "M", "M"), a = c(1, 2, 3)), spec), "stratum 'F'")
})

test_that("stratum-specific supplied cutoffs support sex-specific clinical limits", {
  bm <- list(define_biomarker("hdl", "metabolic", "low", "mg/dL", cutoff = c(Female = 50, Male = 40)),
    define_biomarker("sbp", "cardiovascular", "high", "mmHg", cutoff = 130))
  spec <- define_protocol(bm, "cutoff_count", strata = "sex")
  ref <- fit_reference(spec = spec)
  expect_equal(reference_cutoffs(ref)$lower, c(50, NA, 40, NA))
  d <- data.frame(sex = c("Female", "Male", "Male"), hdl = c(45, 45, 39), sbp = c(120, 135, 120))
  expect_equal(score_load(d, ref)$scores$score, c(1, 1, 1))
  expect_error(define_protocol(bm, "cutoff_count"), "exactly one strata")
  expect_error(define_protocol(list(bm[[2]]), "cutoff_count", strata = "sex"), "stratified cutoffs")
  odd <- define_biomarker("tg", "metabolic", "high", "mg/dL", cutoff = c(Female = 150, Other = 150))
  expect_error(define_protocol(c(bm, list(odd)), "cutoff_count", strata = "sex"), "same stratum levels")
  expect_error(define_biomarker("x", "s", "both", "u", cutoff = list(A = c(2, 1))), "cutoff")
  reduced <- panel_influence(d, ref)$results$hdl$reference
  expect_identical(reduced$spec$strata, character())
  expect_equal(reference_cutoffs(reduced)$upper, 130)
})

test_that("age-conditional quantile references follow the conditional quantile", {
  skip_if_not_installed("quantreg")
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")), age = "age", age_df = 1)
  ref <- fit_reference(age_reference_data(), spec)
  at <- reference_cutoffs(ref, at = data.frame(age = c(30, 50)))
  expect_equal(at$upper, c(131, 151), tolerance = 1e-6)
  expect_equal(nrow(reference_cutoffs(ref)), 5)
  r <- score_load(data.frame(age = c(50, 50, 90), a = c(152, 150, 200)), ref)
  expect_equal(r$scores$score, c(1, 0, NA))
  expect_equal(r$scores$status, c("ok", "ok", "age_outside_reference"))
  wide <- define_protocol(spec$biomarkers, age = "age", age_df = 1, age_extrapolate = TRUE)
  expect_equal(score_load(data.frame(age = 90, a = 200), fit_reference(age_reference_data(), wide))$scores$score, 1)
  expect_error(fit_reference(age_reference_data()[1:15, ], spec), "At least 20")
  expect_error(define_protocol(spec$biomarkers, age = "age", age_df = 11), "age_df")
  expect_error(define_protocol(list(define_biomarker("a", "s", "high", "u", 1)), "cutoff_count", age = "age"), "Age-conditional")
})

test_that("age-conditional z references use the absolute-residual model", {
  age <- rep(20:79, each = 2)
  d <- data.frame(age = age, a = 10 + 0.5 * age + rep(c(-1, 1), times = 60))
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u")), method = "z_mean", age = "age", age_df = 1)
  ref <- fit_reference(d, spec)
  at <- reference_cutoffs(ref, at = data.frame(age = 40))
  expect_equal(at$center, 30, tolerance = 1e-8)
  expect_equal(at$scale, sqrt(pi / 2), tolerance = 1e-8)
  expect_equal(score_load(data.frame(age = 40, a = 30 + 2 * sqrt(pi / 2)), ref)$scores$score, 2, tolerance = 1e-8)
})

test_that("age references validate, compare and plot", {
  skip_if_not_installed("quantreg")
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u"), define_biomarker("b", "t", "low", "u")),
    age = "age", age_df = 2, strata = "sex")
  set.seed(3)
  d <- data.frame(age = runif(400, 20, 80), sex = rep(c("F", "M"), 200))
  d$a <- 100 + 0.5 * d$age + rnorm(400, 0, 5); d$b <- 60 - 0.1 * d$age + rnorm(400, 0, 4)
  ref <- fit_reference(d, spec)
  expect_length(ref$models, 4)
  pooled <- fit_reference(d, define_protocol(spec$biomarkers))
  expect_false(compare_references(list(age = ref, pooled = pooled))$algorithm_equal)
  r <- allostatic_load(d, reference = ref, graphics = "auto")
  expect_match(attr(r$plots, "omitted")$cutoffs, "plot_centiles")
  expect_error(allostaticR:::plot_cutoffs(r), "plot_centiles")
  expect_silent(ggplot2::ggplot_build(plot_centiles(r)))
  expect_silent(ggplot2::ggplot_build(plot_centiles(ref, "a")))
  bad <- ref; bad$models[[1]]$coef$upper <- c(1, NA, 2)
  expect_error(score_load(d, bad), "age model")
  inf <- panel_influence(d, ref)
  expect_equal(nrow(inf$summary), 2)
  boot <- bootstrap_load(d[1:50, ], d, define_protocol(spec$biomarkers, strata = "sex"), B = 5, seed = 1)
  expect_true("stratum" %in% names(boot$cutoff_intervals))
})

test_that("schema 1 passports are upgraded to the current protocol", {
  old_spec <- structure(list(biomarkers = list(structure(list(name = "a", system = "s", direction = "high",
      unit = "u", cutoff = 10, label = "a"), class = "al_biomarker")), method = "cutoff_count",
    risk_quantile = 0.75, missing = "complete", min_coverage = 1, inclusive = TRUE, quantile_type = 7L,
    protocol_id = "custom", protocol_version = "1", sources = character(), medications = character()), class = "al_spec")
  old_ref <- structure(list(spec = old_spec, cutoffs = data.frame(biomarker = "a", label = "a", system = "s",
    direction = "high", unit = "u", lower = NA_real_, upper = 10, n_reference = NA_integer_, stringsAsFactors = FALSE),
    label = "old", n_reference = NA_integer_), class = "al_reference")
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  saveRDS(list(schema_version = 1L, package_version = "0.1.0", protocol = old_ref), f)
  upgraded <- read_protocol(f)
  expect_equal(upgraded$spec$biomarkers[[1]]$transform, "identity")
  expect_equal(score_load(data.frame(a = c(9, 10)), upgraded)$scores$score, c(0, 1))
  saveRDS(list(schema_version = 1L, protocol = old_spec), f)
  expect_s3_class(read_protocol(f), "al_spec")
})

test_that("protocol passports round-trip and reject unsupported schemas", {
  # a reference with fixed cutoffs
  ref <- fit_reference(spec = cutoff_spec())
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  save_protocol(ref, f)
  expect_identical(read_protocol(f), ref)
  # a reference with named biomarkers
  s <- cutoff_spec(); names(s$biomarkers) <- c("marker_a", "marker_b")
  named <- fit_reference(spec = s)
  fn <- tempfile(fileext = ".rds"); on.exit(unlink(fn))
  save_protocol(named, fn)
  expect_identical(read_protocol(fn), named)
  # a reference with age-conditional models
  if (requireNamespace("quantreg", quietly = TRUE)) {
    spec <- define_protocol(list(define_biomarker("a", "s", "high", "u"), define_biomarker("b", "t", "low", "u")),
      age = "age", age_df = 2, strata = "sex")
    set.seed(3)
    d <- data.frame(age = runif(400, 20, 80), sex = rep(c("F", "M"), 200))
    d$a <- 100 + 0.5 * d$age + rnorm(400, 0, 5); d$b <- 60 - 0.1 * d$age + rnorm(400, 0, 4)
    aged <- fit_reference(d, spec)
    fa <- tempfile(fileext = ".rds"); on.exit(unlink(fa))
    save_protocol(aged, fa)
    expect_identical(read_protocol(fa), aged)
  }
})
