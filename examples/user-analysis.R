# Install once before running this script:
# install.packages("allostaticR_0.1.0.tar.gz", repos = NULL, type = "source")
#
# Command-line usage:
# Rscript examples/user-analysis.R data.csv
# Rscript examples/user-analysis.R data.csv reference.csv results
#
# Required columns and units:
# sbp, dbp: mmHg; hdl, triglycerides: mg/dL; hba1c: %; crp: mg/L.
# Optional columns: id (unique participant identifier) and group.
# Measurements must already use the units declared below.

library(allostaticR)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Provide a CSV: Rscript examples/user-analysis.R data.csv")

read_data <- function(path) read.csv(path, na.strings = c("", "NA"), check.names = FALSE)
study_data <- read_data(args[1])
reference_data <- if (length(args) >= 2) read_data(args[2]) else study_data
output_dir <- if (length(args) >= 3) args[3] else "allostaticR_results"
if (file.exists(output_dir)) stop("Output already exists. Choose another directory to preserve its files.")

spec <- define_protocol(
  biomarkers = list(
    define_biomarker("sbp", "cardiovascular", "high", "mmHg", label = "Systolic blood pressure"),
    define_biomarker("dbp", "cardiovascular", "high", "mmHg", label = "Diastolic blood pressure"),
    define_biomarker("hdl", "metabolic", "low", "mg/dL", label = "HDL cholesterol"),
    define_biomarker("triglycerides", "metabolic", "high", "mg/dL", label = "Triglycerides"),
    define_biomarker("hba1c", "metabolic", "high", "%", label = "Glycated hemoglobin"),
    define_biomarker("crp", "inflammatory", "high", "mg/L", label = "C-reactive protein")
  ),
  method = "quantile_count",
  risk_quantile = 0.75,
  missing = "complete",
  inclusive = TRUE
)

# High risk: 75th percentile; HDL: 25th percentile. Equality counts as risk.
# Without reference.csv, thresholds are estimated from the study sample itself.
reference <- fit_reference(reference_data, spec,
  label = if (length(args) >= 2) "External reference" else "Study sample")
result <- allostatic_load(
  study_data, reference = reference,
  id = if ("id" %in% names(study_data)) "id" else NULL,
  group = if ("group" %in% names(study_data)) "group" else NULL,
  graphics = "none", locale = "en"
)

print(summary(result))
print(reference_cutoffs(reference))
if (!dir.create(output_dir, recursive = TRUE)) stop("Could not create the output directory.")
write.csv(as.data.frame(result), file.path(output_dir, "scores.csv"), row.names = FALSE, na = "")
write.csv(score_components(result), file.path(output_dir, "components.csv"), row.names = FALSE, na = "")
write.csv(reference_cutoffs(reference), file.path(output_dir, "cutoffs.csv"), row.names = FALSE, na = "")
# Choose the graph here. All mappings, layers, labels and scales are editable.
plot_data <- plot_data(result, "scores")
if (!"group" %in% names(plot_data)) plot_data$group <- "All participants"
plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = score)) +
  ggplot2::geom_boxplot(width = 0.55, fill = "white", na.rm = TRUE) +
  ggplot2::labs(title = "Allostatic load by study group",
    x = "Study group", y = "Allostatic load score") +
  ggplot2::theme_classic(base_size = 12)
# Examples: add scale_y_continuous(breaks = 0:6), coord_flip(), or theme_bw().
# For a histogram, use ggplot(plot_data, aes(x = score)) + geom_histogram(...).
save_plots(list(group_boxplot = plot), file.path(output_dir, "figures"))
save_protocol(reference, file.path(output_dir, "reference.rds"))
cat("\nResults saved to:", output_dir, "\n")
