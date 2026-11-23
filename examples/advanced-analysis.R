# Run from package root with allostaticR >= 0.1.0 installed.
# Rscript examples/advanced-analysis.R
source("examples/fictitious-analysis.R")
advanced <- file.path(out, "advanced")
dir.create(advanced, recursive = TRUE, showWarnings = FALSE)

spec$protocol_id <- "synthetic-six-marker-panel"
spec$protocol_version <- "1"
spec$sources <- "Demonstration protocol; synthetic quantiles, not clinical cutoffs."
reference <- fit_reference(reference_data, spec, label = "Synthetic independent baseline (n=600)")
result <- allostatic_load(study, reference = reference, id = "id", group = "group", graphics = "auto")

q90 <- spec; q90$risk_quantile <- 0.9; q90$protocol_version <- "2-q90"
partial <- spec; partial$missing <- "proportion"; partial$min_coverage <- 5 / 6
sensitivity <- sensitivity_analysis(study, reference_data,
  list(q75_complete = spec, q90_complete = q90, q75_partial = partial), id = "id")
influence <- panel_influence(study, reference, id = "id")
system_influence <- panel_influence(study, reference, by = "system", id = "id")
bootstrap <- bootstrap_load(study, reference_data, spec, B = 200, seed = 42)
units <- setNames(vapply(biomarkers, `[[`, character(1), "unit"),
  vapply(biomarkers, `[[`, character(1), "name"))
quality <- check_data_quality(study, spec, units = units, id = "id",
  limits = list(sbp = c(60, 240), dbp = c(30, 150), hdl = c(10, 120),
    triglycerides = c(20, 800), hba1c = c(3, 15), crp = c(0, 100)))

visits <- do.call(rbind, lapply(c(0, 6, 12), function(month) {
  d <- study; d$month <- month
  # Artificial improvement in Group B, included only to exercise the API.
  selected <- d$group == "Group B"
  d$sbp[selected] <- d$sbp[selected] - month * 0.5
  d$dbp[selected] <- d$dbp[selected] - month * 0.2
  d$hdl[selected] <- d$hdl[selected] + month * 0.2
  d
}))
trajectory <- score_trajectory(visits, reference, "id", "month", time_unit = "months",
  subjects = head(result$scores$id[result$scores$status == "ok"], 3))
changes <- visit_changes(trajectory)

# Explicit synthetic medication rule; no undocumented general drug adjustment.
set.seed(42)
study$treated_sbp <- runif(nrow(study)) < 0.15
med_spec <- spec; med_spec$medications <- c(sbp = "treated_sbp")
med_result <- score_load(study, fit_reference(reference_data, med_spec), id = "id")
med_comparison <- compare_scores(list(no_rule = result, explicit_rule = med_result))

# Simple mass-concentration conversion: mg/dL -> g/L, supplied factor 0.01.
converted <- convert_units(study, "hdl", "mg/dL", "g/L", factor = 0.01)
conversion_spec <- spec
conversion_spec$biomarkers[[3]]$unit <- "g/L"
converted_reference <- convert_units(reference_data, "hdl", "mg/dL", "g/L", factor = 0.01)
converted_result <- score_load(converted, fit_reference(converted_reference, conversion_spec), id = "id")
stopifnot(identical(converted_result$scores$score, result$scores$score))

passport <- file.path(advanced, "reference-protocol.rds")
save_protocol(reference, passport, overwrite = TRUE)
stopifnot(identical(read_protocol(passport), reference))
compatibility <- compare_references(list(baseline = reference,
  restored = read_protocol(passport), alternative = fit_reference(reference_data, q90)))

write.csv(sensitivity$comparisons, file.path(advanced, "sensitivity.csv"), row.names = FALSE)
write.csv(influence$summary, file.path(advanced, "biomarker-influence.csv"), row.names = FALSE)
write.csv(system_influence$summary, file.path(advanced, "system-influence.csv"), row.names = FALSE)
write.csv(bootstrap$summary, file.path(advanced, "bootstrap-summary.csv"), row.names = FALSE)
write.csv(bootstrap$cutoff_intervals, file.path(advanced, "bootstrap-cutoffs.csv"), row.names = FALSE)
write.csv(changes, file.path(advanced, "longitudinal-changes.csv"), row.names = FALSE)
write.csv(quality$summary, file.path(advanced, "quality.csv"), row.names = FALSE)
write.csv(compatibility, file.path(advanced, "compatibility.csv"), row.names = FALSE)
write.csv(med_comparison, file.path(advanced, "medication-comparison.csv"), row.names = FALSE)
plots <- list(sensitivity = sensitivity$plots$changes,
  biomarker_influence = influence$plots$influence,
  system_influence = system_influence$plots$influence,
  calibration = bootstrap$plots$mean_interval,
  trajectory = trajectory$plots$cohort,
  missingness = quality$plots$missingness,
  cutoffs = plot_cutoffs(result))
save_plots(plots, file.path(advanced, "figures"), height = 150, overwrite = TRUE)
analyses <- list(quality = quality, sensitivity = sensitivity, influence = influence,
  systems = system_influence, bootstrap = bootstrap, trajectory = trajectory)
saveRDS(analyses, file.path(advanced, "analyses.rds"))
render_report(result, file.path(advanced, "report.html"), analyses = analyses,
  locale = "en", overwrite = TRUE)
cat("\nAdvanced example completed: ", advanced, "\n", sep = "")
print(sensitivity$comparisons)
print(bootstrap$summary)
print(compatibility)
