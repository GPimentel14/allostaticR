#' allostaticR: Allostatic Load Estimation and Analysis
#'
#' Start here for the package workflow. Define a biomarker panel, calibrate a
#' reference, calculate scores, and choose how to plot or export the results.
#' All package-generated text is in English. Scoring generates no plots by
#' default; graphics use ordinary editable ggplot2 objects.
#'
#' @section Step 1 - Define the panel:
#' Use [define_biomarker()] to specify a measured column, physiological system,
#' unit, risk direction and optional log transform. Combine biomarkers with
#' [define_protocol()], which fixes the method (`quantile_count`,
#' `cutoff_count` or `z_mean`), the balance (each biomarker or each system
#' weighs equally), the missing-data policy and optional reference `strata`
#' (for example sex) and `age`. HDL normally uses `direction = "low"`; keep
#' its measured values unchanged. Audit inputs with [check_data_quality()].
#'
#' @section Step 2 - Fit a reference and calculate scores:
#' [fit_reference()] estimates fixed thresholds, or reference means and SDs,
#' pooled, within strata, or as smooth functions of age (growth-chart style
#' centiles). Inspect them with [reference_cutoffs()] and [plot_centiles()].
#' Reuse the same reference with [score_load()] to compare cohorts or visits.
#' [allostatic_load()] combines fitting and scoring when supplied a protocol,
#' or scores against a supplied fitted reference.
#'
#' @section Step 3 - Inspect results:
#' `as.data.frame(result)` returns participant scores (use
#' `what = "components"` or `what = "systems"` for biomarker contributions and
#' system scores), and `result$audit` holds calculation metadata.
#' [simulate_biomarkers()] supplies synthetic data, optionally with age and sex.
#'
#' @section Step 4 - Choose and customize a graphic:
#' [plot_data()] returns data for your own `ggplot2::ggplot()` code. [plot_load()]
#' (or `plot(result)`) selects an optional preset; `allostatic_load(..., graphics = "auto")`
#' builds the full preset collection (including cutoff panels via
#' `plot_load(result, "cutoffs")`). [plot_centiles()] shows age-specific curves
#' and [plot_risk_map()] individual contributions. [theme_allostatic()] provides
#' an optional theme.
#'
#' @section Step 5 - Explain differences and patterns:
#' [decompose_gap()] splits a difference in mean load between two groups
#' exactly into biomarker or system contributions. [risk_patterns()] shows
#' which combinations of biomarkers make up each load, their co-occurrence and
#' the redundancy of the panel.
#'
#' @section Step 6 - Quantify uncertainty and robustness:
#' [measurement_uncertainty()] propagates analytical and within-person
#' biological variation to each score; [visit_changes()] can classify changes
#' with reference change values. [bootstrap_load()] measures
#' reference-calibration uncertainty. [sensitivity_analysis()],
#' [panel_influence()], [compare_scores()] and the multiverse tools
#' [multiverse_grid()] and [multiverse_analysis()] assess scoring choices.
#'
#' @section Step 7 - Handle missing data explicitly:
#' Complete panels are the default. [check_data_quality()] audits inputs and
#' [convert_units()] applies an explicit factor and offset.
#' [impute_biomarkers()] performs explicit multiple imputation (random
#' forests, missForest, predictive mean matching and others),
#' [score_imputed()] pools scores with Rubin's rules and [check_imputation()]
#' compares methods on masked values.
#'
#' @section Step 8 - Follow-up, save and reproduce:
#' [score_trajectory()] uses a fixed reference across visits;
#' [visit_changes()] calculates paired endpoint changes. [save_protocol()] and
#' [read_protocol()] persist validated protocol passports.
#' [compare_references()] checks algorithm and threshold compatibility.
#' [model_trajectory()] fits mixed models to trajectories.
#'
#' @section Step 9 - Surveys, validation and study design:
#' Protocol `weights` and [survey_summary()] support complex surveys;
#' [validate_scores()] compares algorithms against outcomes with
#' cross-validation; [protocol_template()] rebuilds published panels;
#' `medication_rule = "add"` applies treatment corrections;
#' [harmonize_cohorts()] aligns cohorts; [detect_outliers()] screens values;
#' [tipping_point()] tests sensitivity to data missing not at random.
#'
#' @section Where to find help:
#' Use `help("allostaticR")` for this overview, `help("allostatic_load")`
#' for scoring, and `help("plot_data")` for the customizable plotting example.
#' `vignette("introduction", package = "allostaticR")` provides a worked
#' tutorial and `vignette("advanced-methods", package = "allostaticR")` covers
#' normative references, decomposition, patterns, uncertainty, multiverse
#' analysis and imputation.
#'
#' @examples
#' source(system.file("examples", "panel.R", package = "allostaticR"))
#' reference <- fit_reference(data, spec, label = "Study baseline")
#' result <- score_load(data, reference, id = "id", group = "group")
#' summary(result)
#' scores <- plot_data(result)
#' # See ?plot_data for a complete, editable ggplot2 boxplot example.
#'
#' @keywords internal
#' @importFrom ggplot2 ggplot
"_PACKAGE"

