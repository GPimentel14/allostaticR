# Canonical three-biomarker panel used by the package examples, vignettes
# and README. Sourcing this file defines data and spec.
data <- simulate_biomarkers(100)
spec <- define_protocol(list(
  define_biomarker("sbp", "cardiovascular", "high", "mmHg"),
  define_biomarker("hdl", "metabolic", "low", "mg/dL"),
  define_biomarker("crp", "inflammatory", "high", "mg/L")
))
