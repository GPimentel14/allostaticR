# Shared fixtures for the test suite.
# Every helper here is deterministic: the same call always returns the same
# object, so tests never depend on the order in which they run.

# Two biomarkers with fixed cutoffs; the workhorse protocol of the suite.
cutoff_spec <- function(missing = "complete", min_coverage = 1, inclusive = TRUE) {
  define_protocol(list(define_biomarker("a", "s1", "high", "u", cutoff = 10),
               define_biomarker("b", "s2", "low", "u", cutoff = 5)),
          method = "cutoff_count", missing = missing,
          min_coverage = min_coverage, inclusive = inclusive)
}

# The canonical three-biomarker panel (sbp/hdl/crp) used across the package.
panel_spec <- function() {
  define_protocol(list(define_biomarker("sbp", "cardiovascular", "high", "mmHg"),
    define_biomarker("hdl", "metabolic", "low", "mg/dL"),
    define_biomarker("crp", "inflammatory", "high", "mg/L")))
}

# The same panel under the name used by the multiverse domain.
multiverse_spec <- function() panel_spec()

# Synthetic participants for the canonical panel.
panel_data <- function(n = 100L, seed = 1L) simulate_biomarkers(n, seed = seed)

# Participants with missing values in two of the three biomarkers.
imputation_data <- function(n = 120) {
  set.seed(11)
  d <- data.frame(id = seq_len(n), sex = rep(c("F", "M"), n / 2), a = stats::rnorm(n, 10, 2))
  d$b <- 2 * d$a + stats::rnorm(n, 0, 2); d$c <- exp(stats::rnorm(n))
  d$b[sample.int(n, 20)] <- NA; d$c[sample.int(n, 15)] <- NA
  d
}

# Protocol of the imputation fixtures: two markers of one system and a
# log-transformed marker of another one.
imputation_spec <- function(...) define_protocol(list(define_biomarker("a", "s1", "high", "u"),
  define_biomarker("b", "s1", "high", "u"), define_biomarker("c", "s2", "high", "u", transform = "log")), ...)

# Sample used to fit age-conditional references.
age_reference_data <- function() {
  age <- rep(20:79, each = 5)
  data.frame(age = age, a = 100 + age + rep(c(-2, -1, 0, 1, 2), times = 60))
}

# Coefficients of variation for the measurement-uncertainty model.
cv_table <- function(cva = 0.02, cvi = 0.1, markers = c("a", "b")) {
  data.frame(biomarker = markers, cv_analytical = cva, cv_within = cvi, source = "test values")
}

# Exact distribution of a weighted sum of independent Bernoulli variables.
brute_force <- function(p, weights) {
  combos <- as.matrix(expand.grid(rep(list(0:1), length(p))))
  prob <- apply(combos, 1, function(k) prod(ifelse(k == 1, p, 1 - p)))
  value <- round(drop(combos %*% weights), 10)
  out <- tapply(prob, value, sum)
  list(value = as.numeric(names(out)), prob = as.vector(out))
}

# A scored result, ready to be reused by several domains.
scored_result <- function(spec = cutoff_spec(), data = NULL) {
  if (is.null(data)) data <- data.frame(a = c(0, 10, 10, NA_real_), b = c(5, 6, 4, 4))
  allostatic_load(data, spec)
}
