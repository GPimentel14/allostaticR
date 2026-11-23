# Scoring: score_load(), allostatic_load(), as.data.frame.al_result(),
# the medication rules and simulate_biomarkers().

test_that("additive medication effects correct measured values", {
  b <- list(define_biomarker("sbp", "cv", "high", "mmHg", cutoff = 140))
  spec <- define_protocol(b, "cutoff_count", medications = c(sbp = "treated"), medication_rule = "add",
    medication_effects = c(sbp = 10))
  r <- score_load(data.frame(sbp = c(135, 135, 125), treated = c(FALSE, TRUE, TRUE)), fit_reference(spec = spec))
  expect_equal(r$scores$score, c(0, 1, 0))
  expect_equal(as.data.frame(r, what = "components")$value_adjusted, c(135, 145, 135))
  expect_equal(as.data.frame(r, what = "components")$value, c(135, 135, 125))
  zspec <- define_protocol(list(define_biomarker("sbp", "cv", "high", "mmHg")), method = "z_mean",
    medications = c(sbp = "treated"), medication_rule = "add", medication_effects = c(sbp = 10))
  ref <- fit_reference(data.frame(sbp = c(120, 130, 120), treated = c(FALSE, FALSE, TRUE)), zspec)
  expect_equal(reference_cutoffs(ref)$center, mean(c(120, 130, 130)))
  expect_error(define_protocol(b, "cutoff_count", medications = c(sbp = "t"), medication_rule = "add"), "medication_effects")
  expect_error(define_protocol(list(define_biomarker("sbp", "cv", "high", "u")), method = "z_mean", medications = c(sbp = "t")), "medication_rule")
})

test_that("medication rules are explicit and never impute missing measurements", {
  s <- cutoff_spec(); s$medications <- c(a = "treated")
  d <- data.frame(a = c(0, 10, NA), b = c(6, 6, 6), treated = c(TRUE, TRUE, TRUE))
  r <- score_load(d, fit_reference(spec = s))
  expect_equal(r$scores$score, c(1, 1, NA))
  expect_equal(r$components$medication_adjusted[1:3], c(TRUE, FALSE, FALSE))
  d$treated[1] <- NA
  expect_error(score_load(d, fit_reference(spec = s)), "complete logical")
})

test_that("incompatible reference metadata and unused cuts cannot change scoring silently", {
  ref <- fit_reference(spec = cutoff_spec()); ref$cutoffs$lower[1] <- 0
  expect_error(score_load(data.frame(a = 0, b = 6), ref), "contradicts")
  ref <- fit_reference(spec = cutoff_spec()); ref$cutoffs$unit[1] <- "other"
  expect_error(score_load(data.frame(a = 0, b = 6), ref), "metadata")
  spec <- cutoff_spec(); spec$medications <- c(a = "treated")
  q <- check_data_quality(data.frame(a = 0, b = 6), spec)
  expect_true("invalid_medication_column:treated" %in% q$issues)
})

test_that("cutoff directions, equality, identifiers and missingness are explicit", {
  d <- data.frame(id = c("z", "x", "y", "w"), a = c(10, 9, 11, NA), b = c(5, 6, 4, 4))
  ref <- fit_reference(spec = cutoff_spec())
  r <- score_load(d, ref, id = "id")
  expect_equal(r$scores$score, c(2, 0, 2, NA))
  expect_equal(r$scores$id, d$id)
  expect_equal(r$scores$n_observed, c(2, 2, 2, 1))
  expect_equal(r$scores$status[4], "insufficient_coverage")
  expect_equal(as.data.frame(r, what = "components")$risk, c(1L, 0L, 1L, NA, 1L, 0L, 1L, 1L))
  strict <- score_load(d, fit_reference(spec = cutoff_spec(inclusive = FALSE)))
  expect_equal(strict$scores$score, c(0, 0, 2, NA))
  prop <- score_load(d, fit_reference(spec = cutoff_spec("proportion", 0.5)))
  expect_equal(prop$scores$score, c(1, 0, 1, 1))
  expect_equal(prop$audit$scale, "proportion")
})

test_that("reference quantiles remain fixed when applying new cohorts", {
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u"),
                       define_biomarker("b", "s", "low", "u")))
  ref <- fit_reference(data.frame(a = 1:4, b = 1:4), spec, "external")
  before <- ref
  expect_equal(reference_cutoffs(ref)$upper[1], 3.25)
  expect_equal(reference_cutoffs(ref)$lower[2], 1.75)
  r <- score_load(data.frame(a = c(100, 2), b = c(100, 1)), ref)
  expect_equal(r$scores$score, c(1, 1))
  expect_identical(ref, before)
  expect_identical(r$reference, before)
  expect_equal(as.data.frame(r, what = "systems")$risk_proportion, c(0.5, 0.5))
})

test_that("bilateral thresholds and ties obey protocol", {
  spec <- define_protocol(list(define_biomarker("a", "s", "both", "u", c(2, 8))), "cutoff_count")
  r <- score_load(data.frame(a = c(1, 2, 5, 8, 9, NA)), fit_reference(spec = spec))
  expect_equal(r$scores$score, c(1, 1, 0, 1, 1, NA))
  qspec <- define_protocol(list(define_biomarker("a", "s", "high", "u")))
  ties <- allostatic_load(data.frame(a = rep(1, 4)), qspec, graphics = "none")
  expect_equal(ties$scores$score, rep(1, 4))
})

test_that("all missing panels cannot yield zero risk", {
  d <- data.frame(a = NA_real_, b = NA_real_)
  r <- score_load(d, fit_reference(spec = cutoff_spec("proportion", 0.5)))
  expect_true(is.na(r$scores$score))
  expect_true(all(is.na(r$systems$risk_proportion)))
})

test_that("scoring rejects invalid inputs", {
  ref <- fit_reference(spec = cutoff_spec())
  expect_error(score_load(data.frame(a = 1), ref), "Missing biomarker")
  expect_error(score_load(data.frame(a = "1", b = 2), ref), "numeric")
  expect_error(score_load(data.frame(a = Inf, b = 2), ref), "non-finite")
  expect_error(score_load(data.frame(a = NaN, b = 2), ref), "non-finite")
  expect_error(score_load(data.frame(a = 1:2, b = 2:3, id = c(1, 1)), ref, id = "id"), "unique")
  expect_error(allostatic_load(data.frame(a = 1, b = 2), cutoff_spec(), ref), "not both")
})

test_that("z_mean uses fixed reference means and SDs with risk orientation", {
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"),
    define_biomarker("b", "s2", "low", "u")), method = "z_mean")
  ref <- fit_reference(data.frame(a = 1:5, b = c(10, 20, 30, 40, 50)), spec)
  expect_equal(reference_cutoffs(ref)$center, c(3, 30))
  expect_equal(reference_cutoffs(ref)$scale, c(sqrt(2.5), sqrt(250)))
  r <- score_load(data.frame(a = c(3 + sqrt(2.5), 3), b = c(30, 30 - 2 * sqrt(250))), ref)
  expect_equal(r$scores$score, c(0.5, 1))
  expect_equal(r$audit$scale, "z")
  expect_equal(as.data.frame(r, what = "components")$z, c(1, 0, 0, 2))
  expect_true(all(is.na(as.data.frame(r, what = "components")$risk)))
  expect_equal(as.data.frame(r, what = "systems")$score, c(1, 0, 0, 2))
})

test_that("log transforms and caps apply to z-scores only", {
  spec <- define_protocol(list(define_biomarker("a", "s", "high", "u", transform = "log")), method = "z_mean")
  ref <- fit_reference(data.frame(a = exp(0:2)), spec)
  expect_equal(score_load(data.frame(a = exp(3)), ref)$scores$score, 2)
  capped <- fit_reference(data.frame(a = exp(0:2)), define_protocol(spec$biomarkers, method = "z_mean", z_cap = 1.5))
  expect_equal(score_load(data.frame(a = exp(c(3, -3))), capped)$scores$score, c(1.5, -1.5))
  expect_error(score_load(data.frame(a = 0), ref), "positive")
  q <- define_protocol(list(define_biomarker("a", "s", "high", "u", transform = "log")))
  expect_equal(reference_cutoffs(fit_reference(data.frame(a = 1:4), q))$upper, 3.25)
})

test_that("system balance gives every system equal weight", {
  bm <- list(define_biomarker("a", "s1", "high", "u", 10), define_biomarker("b", "s1", "high", "u", 10),
    define_biomarker("c", "s2", "high", "u", 10))
  d <- data.frame(a = c(11, 11), b = c(0, 11), c = c(11, 0))
  balanced <- score_load(d, fit_reference(spec = define_protocol(bm, "cutoff_count", balance = "system")))
  expect_equal(balanced$scores$score, c(0.75, 0.5))
  expect_equal(balanced$audit$scale, "system_proportion")
  expect_equal(score_load(d, fit_reference(spec = define_protocol(bm, "cutoff_count")))$scores$score, c(2, 2))
  partial <- define_protocol(bm, "cutoff_count", balance = "system", missing = "proportion", min_coverage = 0.3)
  r <- score_load(data.frame(a = c(NA, 11), b = c(11, NA), c = c(0, NA)), fit_reference(spec = partial))
  expect_equal(r$scores$score, c(0.5, NA))
  expect_equal(r$scores$status, c("ok", "insufficient_coverage"))
  zspec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s1", "high", "u"),
    define_biomarker("c", "s2", "high", "u")), method = "z_mean", balance = "system")
  ref <- fit_reference(data.frame(a = c(-1, 1), b = c(-1, 1), c = c(-1, 1)), zspec)
  s <- sqrt(2)
  z <- score_load(data.frame(a = s, b = 3 * s, c = 0), ref)
  expect_equal(z$scores$score, 1)
  expect_equal(z$audit$scale, "system_z")
})

test_that("bilateral quantile tails can be split", {
  b <- define_biomarker("a", "s", "both", "u")
  full <- reference_cutoffs(fit_reference(data.frame(a = 1:100), define_protocol(list(b))))
  split <- reference_cutoffs(fit_reference(data.frame(a = 1:100), define_protocol(list(b), bilateral = "split")))
  expect_equal(c(full$lower, full$upper), unname(stats::quantile(1:100, c(0.25, 0.75))))
  expect_equal(c(split$lower, split$upper), unname(stats::quantile(1:100, c(0.125, 0.875))))
  flagged <- score_load(data.frame(a = 1:100), fit_reference(data.frame(a = 1:100),
    define_protocol(list(b), bilateral = "split")))
  expect_equal(mean(flagged$scores$score), 0.26)
})

test_that("summaries, references and versions print and record metadata", {
  r <- allostatic_load(data.frame(a = c(0, 10), b = c(5, 6)), cutoff_spec())
  s <- summary(r)
  expect_s3_class(s, "summary.al_result")
  expect_output(print(s), "Allostatic load summary")
  expect_output(print(r$reference), "Allostatic load reference")
  expect_equal(r$audit$package_version, as.character(utils::packageVersion("allostaticR")))
  expect_identical(as.data.frame(r, what = "components")$contribution, as.numeric(as.data.frame(r, what = "components")$risk))
})

test_that("simulation is reproducible without changing caller RNG", {
  set.seed(75); before <- .Random.seed
  expect_identical(simulate_biomarkers(5, 1), simulate_biomarkers(5, 1))
  expect_identical(.Random.seed, before)
  expect_error(simulate_biomarkers(0), "positive integer")
})

test_that("default simulation is unchanged and covariates are optional", {
  d <- simulate_biomarkers(3, seed = 1)
  expect_named(d, c("id", "group", "sbp", "hdl", "crp"))
  set.seed(1); expect_equal(d$sbp, round(stats::rnorm(3, 125, 15), 1))
  x <- simulate_biomarkers(50, seed = 2, covariates = TRUE)
  expect_true(all(c("age", "sex") %in% names(x)))
  expect_true(all(x$age >= 20 & x$age <= 85))
})
