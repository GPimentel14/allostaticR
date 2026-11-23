# Run from the package root after installing allostaticR:
# Rscript examples/fictitious-analysis.R [output_directory]
library(allostaticR)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "artifacts/demo"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261003)

# Independent biomarker draws for illustration, not a physiological model.
simulate_panel <- function(n) {
  data.frame(
    id = seq_len(n),
    sbp = round(rnorm(n, 122, 14), 1),
    dbp = round(rnorm(n, 77, 9), 1),
    hdl = round(pmax(15, rnorm(n, 55, 12)), 1),
    triglycerides = round(rlnorm(n, log(120), 0.35), 1),
    hba1c = round(pmax(3.5, rnorm(n, 5.4, 0.45)), 2),
    crp = round(rlnorm(n, log(1.2), 0.7), 2)
  )
}
reference_data <- simulate_panel(600)
study <- simulate_panel(400)
study$group <- rep(c("Group A", "Group B"), each = 200)
b <- study$group == "Group B"
# Deliberately introduce differences; these are not estimates of real effects.
study$sbp[b] <- study$sbp[b] + 12
study$dbp[b] <- study$dbp[b] + 5
study$hdl[b] <- study$hdl[b] - 6
study$triglycerides[b] <- round(study$triglycerides[b] * 1.3, 1)
study$hba1c[b] <- study$hba1c[b] + 0.4
study$crp[b] <- round(study$crp[b] * 1.5, 2)
biomarkers <- list(
  define_biomarker("sbp", "cardiovascular", "high", "mmHg", label = "Systolic blood pressure"),
  define_biomarker("dbp", "cardiovascular", "high", "mmHg", label = "Diastolic blood pressure"),
  define_biomarker("hdl", "metabolic", "low", "mg/dL", label = "HDL cholesterol"),
  define_biomarker("triglycerides", "metabolic", "high", "mg/dL", label = "Triglycerides"),
  define_biomarker("hba1c", "metabolic", "high", "%", label = "Glycated hemoglobin"),
  define_biomarker("crp", "inflammatory", "high", "mg/L", label = "C-reactive protein")
)
# Exactly 5% missing per biomarker, selected independently of group and value.
for (marker in biomarkers) study[sample.int(nrow(study), 20), marker$name] <- NA_real_
spec <- define_protocol(biomarkers, method = "quantile_count", missing = "complete")
reference <- fit_reference(reference_data, spec, label = "Independent synthetic reference (n=600)")
result <- allostatic_load(study, reference = reference, id = "id", group = "group", graphics = "auto")

# Independently reconstruct counts to verify the package's output.
manual <- vapply(seq_along(biomarkers), function(j) {
  x <- study[[biomarkers[[j]]$name]]
  if (biomarkers[[j]]$direction == "low") as.integer(x <= reference$cutoffs$lower[j])
  else as.integer(x >= reference$cutoffs$upper[j])
}, integer(nrow(study)))
stopifnot(identical(as.numeric(rowSums(manual)), result$scores$score))

group_summary <- do.call(rbind, lapply(split(result$scores, result$scores$group), function(d) {
  x <- d$score[!is.na(d$score)]
  data.frame(group = d$group[1], n = nrow(d), n_valid = length(x),
    n_missing = sum(is.na(d$score)), mean = mean(x), sd = sd(x),
    median = median(x), q25 = unname(quantile(x, 0.25)), q75 = unname(quantile(x, 0.75)))
}))
rownames(group_summary) <- NULL
components <- score_components(result)
biomarker_summary <- do.call(rbind, lapply(split(components, components$biomarker), function(d) {
  data.frame(biomarker = d$biomarker[1], n_observed = sum(!is.na(d$risk)),
    n_missing = sum(is.na(d$risk)), proportion_at_risk = mean(d$risk, na.rm = TRUE))
}))
# Alternative missing-data policy: fraction at risk, requiring at least 5/6 markers.
proportion_spec <- define_protocol(biomarkers, missing = "proportion", min_coverage = 5 / 6)
proportion_reference <- fit_reference(reference_data, proportion_spec, label = reference$label)
proportion_result <- allostatic_load(study, reference = proportion_reference,
  id = "id", group = "group", graphics = "none")
stopifnot(identical(reference_cutoffs(reference), reference_cutoffs(proportion_reference)))

write.csv(study, file.path(out, "fictitious-study.csv"), row.names = FALSE, na = "")
write.csv(reference_data, file.path(out, "fictitious-reference.csv"), row.names = FALSE)
write.csv(reference_cutoffs(reference), file.path(out, "cutoffs.csv"), row.names = FALSE)
write.csv(as.data.frame(result), file.path(out, "scores.csv"), row.names = FALSE, na = "")
write.csv(as.data.frame(proportion_result), file.path(out, "proportion-scores.csv"), row.names = FALSE, na = "")
write.csv(components, file.path(out, "components.csv"), row.names = FALSE, na = "")
write.csv(group_summary, file.path(out, "group-summary.csv"), row.names = FALSE)
write.csv(biomarker_summary, file.path(out, "biomarker-summary.csv"), row.names = FALSE)
save_plots(result$plots, file.path(out, "figures"), height = 140, overwrite = TRUE)
saveRDS(list(result = result, proportion_result = proportion_result), file.path(out, "analysis.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

report <- c(
  paste("# Demonstration analysis with allostaticR", packageVersion("allostaticR")), "",
  "Entirely synthetic data. Seed: 20261003. These data do not represent a real population.", "",
  "## Design", "",
  "400 participants, two groups of 200 and six biomarkers. Independent reference of 600 records.",
  "Thresholds: 75th percentile for high risk and 25th percentile for HDL (low risk); equality included.",
  "Artificial differences were introduced in all Group B biomarkers, with 5% missingness per biomarker.",
  "Biomarkers were simulated independently; physiological correlations are not reproduced.", "",
  "## Primary result", "",
  sprintf("Complete panel: %d of 400 participants scored; %d unscored due to missing data.",
    sum(!is.na(result$scores$score)), sum(is.na(result$scores$score))),
  sprintf("Overall mean load: %.2f out of 6; median: %.1f.", mean(result$scores$score, na.rm = TRUE), median(result$scores$score, na.rm = TRUE)),
  sprintf("%s: n scored=%d, mean=%.2f, median=%.1f, IQR=%.1f-%.1f.",
    group_summary$group, group_summary$n_valid, group_summary$mean, group_summary$median,
    group_summary$q25, group_summary$q75),
  sprintf("Descriptive difference B minus A: %.2f points among participants with complete panels.",
    group_summary$mean[2] - group_summary$mean[1]), "",
  "## Alternative missing-data policy", "",
  sprintf("With at least 5 of 6 biomarkers: %d scored; mean at-risk proportion=%.3f.",
    sum(!is.na(proportion_result$scores$score)), mean(proportion_result$scores$score, na.rm = TRUE)),
  "This proportion uses a 0-1 scale and should not be compared directly with the 0-6 count.", "",
  "## Interpretation and verification", "",
  "The group difference reflects changes imposed during simulation; it does not establish causality or clinical validity.",
  "No clinical categories were assigned and no significance tests were performed.",
  "Counts were reconstructed independently and match the package output, including NA values.",
  "The proportion alternative uses exactly the same reference thresholds.", "",
  "## Graphics", "",
  "![Distribution](figures/distribution.png)", "",
  "![Groups](figures/groups.png)", "",
  "![Biomarkers](figures/biomarkers.png)", "",
  "![Coverage](figures/coverage.png)", "",
  "## Reproduction", "", "From the package root, with allostaticR installed:", "",
  "```sh", "Rscript examples/fictitious-analysis.R", "```",
  "", "CSV files, RDS objects, PDF/PNG graphics and session-info.txt are saved in artifacts/demo."
)
writeLines(report, file.path(out, "REPORT.md"))
print(group_summary)
cat("\n", paste(report[grep("Complete panel:|Overall mean|Descriptive difference|With at least", report)], collapse = "\n"), "\n", sep = "")
