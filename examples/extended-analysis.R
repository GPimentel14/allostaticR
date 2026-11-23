# Run from the package root with the current allostaticR installed (and the
# optional quantreg, patchwork, mice and randomForest packages):
# Rscript examples/extended-analysis.R
# Reproduces the README figures (man/figures) and writes tables, figures and an
# HTML report to artifacts/demo/extended. All data are synthetic.
source("examples/fictitious-analysis.R")
extended <- file.path(out, "extended")
dir.create(extended, recursive = TRUE, showWarnings = FALSE)
readme_figures <- file.path("man", "figures")
dir.create(readme_figures, recursive = TRUE, showWarnings = FALSE)
readme_png <- function(plot, name, width, height) ggplot2::ggsave(file.path(readme_figures, paste0(name, ".png")),
  plot, width = width, height = height, units = "mm", dpi = 110, bg = "white")

# 1. Normative references by sex and age (simulated data with covariates).
normative_reference_data <- simulate_biomarkers(3000, seed = 11, covariates = TRUE)
normative_study <- simulate_biomarkers(600, seed = 12, covariates = TRUE)
normative_markers <- list(
  define_biomarker("sbp", "cardiovascular", "high", "mmHg", label = "Systolic blood pressure"),
  define_biomarker("hdl", "metabolic", "low", "mg/dL", label = "HDL cholesterol"),
  define_biomarker("crp", "inflammatory", "high", "mg/L", label = "C-reactive protein", transform = "log"))
normative_spec <- define_protocol(normative_markers, strata = "sex", age = "age",
  protocol_id = "synthetic-normative-panel", sources = "Synthetic demonstration; not clinical centiles.")
normative_reference <- fit_reference(normative_reference_data, normative_spec, label = "Synthetic reference (n=3000)")
normative_result <- score_load(normative_study, normative_reference, id = "id", group = "group")
centiles <- plot_centiles(normative_result)
pooled_result <- score_load(normative_study, fit_reference(normative_reference_data,
  define_protocol(normative_markers)), id = "id", group = "group")
age_band <- cut(normative_study$age, c(20, 40, 60, 86), right = FALSE, labels = c("20-39", "40-59", "60-85"))
flag_rates <- data.frame(age_band = levels(age_band),
  pooled_mean_load = tapply(pooled_result$scores$score, age_band, mean),
  normative_mean_load = tapply(normative_result$scores$score, age_band, mean, na.rm = TRUE))
write.csv(flag_rates, file.path(extended, "load-by-age-band.csv"), row.names = FALSE)

# 2. Decomposition of the group difference (six-marker demonstration study).
gap <- decompose_gap(result, B = 1000, seed = 1)
gap_by_system <- decompose_gap(result, by = "system", B = 1000, seed = 1)

# 3. Risk patterns and risk map.
patterns <- risk_patterns(result)
risk_map <- plot_risk_map(result)

# 4. Measurement and biological variation. Illustrative CVs only.
variation <- data.frame(biomarker = vapply(biomarkers, `[[`, character(1), "name"),
  cv_analytical = c(0.02, 0.02, 0.02, 0.03, 0.02, 0.05),
  cv_within = c(0.08, 0.08, 0.07, 0.20, 0.02, 0.40),
  source = "Illustrative values for synthetic data; use sourced CVs (e.g. EFLM database) in real analyses")
uncertainty <- measurement_uncertainty(result, variation)

# 5. Multiverse with a synthetic outcome.
filled <- study
for (b in vapply(biomarkers, `[[`, character(1), "name")) filled[[b]][is.na(filled[[b]])] <- median(filled[[b]], na.rm = TRUE)
set.seed(5)
study$walking_speed <- 1.3 - 0.004 * (filled$sbp - 120) - 0.05 * (filled$hba1c - 5.4) - 0.02 * log(filled$crp) + rnorm(nrow(study), 0, 0.15)
grid <- multiverse_grid(spec, risk_quantile = c(0.75, 0.8, 0.9), missing = c("complete", "proportion"), min_coverage = 5 / 6)
multiverse <- multiverse_analysis(study, reference_data, grid, id = "id", group = "group", outcome = "walking_speed")

# 6. Explicit multiple imputation of the 5% missing measurements.
imputation <- impute_biomarkers(study, spec, method = "rf", m = 5, predictors = "group", seed = 7)
pooled <- score_imputed(imputation, reference, id = "id", group = "group")
imputation_check <- check_imputation(study, spec, reference, methods = c("rf", "missforest", "pmm", "cart", "norm"),
  predictors = "group", prop = 0.15, m = 5, seed = 3)

write.csv(gap$contributions, file.path(extended, "decomposition.csv"), row.names = FALSE)
write.csv(gap_by_system$contributions, file.path(extended, "decomposition-systems.csv"), row.names = FALSE)
write.csv(patterns$patterns, file.path(extended, "risk-patterns.csv"), row.names = FALSE)
write.csv(patterns$cooccurrence, file.path(extended, "cooccurrence.csv"), row.names = FALSE)
write.csv(uncertainty$biomarkers, file.path(extended, "uncertainty-biomarkers.csv"), row.names = FALSE)
write.csv(uncertainty$summary, file.path(extended, "uncertainty-summary.csv"), row.names = FALSE)
write.csv(multiverse$choices, file.path(extended, "multiverse-specifications.csv"), row.names = FALSE)
write.csv(multiverse$effects, file.path(extended, "multiverse-effects.csv"), row.names = FALSE)
write.csv(pooled$pooled, file.path(extended, "imputation-pooled.csv"), row.names = FALSE)
write.csv(imputation_check$summary, file.path(extended, "imputation-check.csv"), row.names = FALSE)

save_plots(list(centiles = centiles, waterfall = gap$plots$waterfall, waterfall_systems = gap_by_system$plots$waterfall,
    risk_map = risk_map, patterns = patterns$plots$patterns, cooccurrence = patterns$plots$cooccurrence,
    uncertainty = uncertainty$plots$reclassification, borderline = uncertainty$plots$borderline,
    stability = multiverse$plots$stability, specification_curve = multiverse$plots$specification_curve,
    imputed_values = imputation$plots$distributions, imputation_check = imputation_check$plots$flags,
    pooled_means = pooled$plots$pooled_means),
  file.path(extended, "figures"), width = 230, height = 160, overwrite = TRUE)
readme_png(centiles, "centiles", 210, 200)
readme_png(gap$plots$waterfall, "decomposition", 220, 125)
readme_png(risk_map, "risk-map", 260, 145)
readme_png(patterns$plots$patterns, "patterns", 220, 150)
readme_png(uncertainty$plots$reclassification, "uncertainty", 170, 140)
readme_png(multiverse$plots$specification_curve, "multiverse", 210, 170)
readme_png(imputation_check$plots$flags, "imputation", 200, 135)
render_report(result, file.path(extended, "report.html"), overwrite = TRUE,
  analyses = list(decomposition = gap, systems = gap_by_system, patterns = patterns, uncertainty = uncertainty,
    multiverse = multiverse, imputation = imputation, imputation_check = imputation_check, pooled = pooled),
  plots = list(risk_map = risk_map))
cat("\nExtended example completed: ", extended, "\n", sep = "")
print(gap); print(patterns); print(uncertainty); print(multiverse); print(pooled); print(imputation_check)
