# Protocol definition and validation: define_biomarker(), define_protocol(),
# allostaticR:::validate_protocol(), allostaticR:::validate_data(), the cutoff rules and the templates.

test_that("published templates rebuild documented panels", {
  t <- protocol_template()
  expect_equal(sum(t$template == "seeman1997"), 10)
  spec <- protocol_template("seeman1997", columns = c(dhea_s = "dheas"), drop = "urinary_epinephrine")
  names_all <- vapply(spec$biomarkers, `[[`, character(1), "name")
  expect_length(names_all, 9); expect_true("dheas" %in% names_all)
  expect_match(spec$sources[1], "Seeman")
  expect_equal(spec$biomarkers[[which(names_all == "hdl")]]$direction, "low")
  expect_error(protocol_template("unknown"), "Unknown template")
})

test_that("z_mean rejects undefined continuous rules", {
  both <- define_biomarker("a", "s", "both", "u")
  expect_error(define_protocol(list(both), method = "z_mean"), "high and low")
  high <- define_biomarker("a", "s", "high", "u")
  expect_error(define_protocol(list(high), method = "z_mean", medications = c(a = "m")), "count method")
  expect_error(define_protocol(list(high), z_cap = 0), "z_cap")
  expect_error(fit_reference(data.frame(a = c(2, 2)), define_protocol(list(high), method = "z_mean")), "SD")
})

test_that("protocol definition rejects invalid inputs", {
  expect_error(define_biomarker("a", "s", "both", "u", c(2, 1)), "cutoff")
  expect_error(define_protocol(list()), "non-empty")
  b <- define_biomarker("a", "s", "high", "u")
  expect_error(define_protocol(list(b, b)), "unique")
  expect_error(define_protocol(list(b), "cutoff_count"), "cutoff")
  expect_error(define_protocol(list(b), risk_quantile = NA_real_), "risk_quantile")
  expect_error(fit_reference(data.frame(a = NA_real_), define_protocol(list(b))), "two observed")
})

test_that("named biomarker lists are supported by the protocol", {
  s <- cutoff_spec(); names(s$biomarkers) <- c("marker_a", "marker_b")
  ref <- fit_reference(spec = s)
  expect_equal(score_load(data.frame(a = 10, b = 5), ref)$scores$score, 2)
})

test_that("protocol and data validation report the problems they find", {
  spec <- cutoff_spec()
  expect_true(allostaticR:::validate_protocol(spec))
  expect_true(allostaticR:::validate_data(data.frame(a = 1, b = 2), spec))
  expect_error(allostaticR:::validate_protocol(list()), "define_protocol")
  expect_error(allostaticR:::validate_data(data.frame(a = 1), spec), "Missing biomarker column")
  expect_error(allostaticR:::validate_data(data.frame(a = "1", b = 2), spec), "numeric vector")
  expect_error(allostaticR:::validate_data(data.frame(a = Inf, b = 2), spec), "non-finite")
  wrong_units <- data.frame(a = 1, b = 2)
  attr(wrong_units, "al_units") <- c(a = "other")
  expect_error(allostaticR:::validate_data(wrong_units, spec), "Declared unit mismatch")
})
