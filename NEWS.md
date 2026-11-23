# allostaticR 0.1.0

## Maintenance

* Internal consolidation: a single rule for risk flags, shared by scoring and
  by the imputation check; one cluster-resampling routine, shared by the
  bootstrap and the survey summary; and shared helpers for biomarker
  metadata, captions, palette, table printing and figure sampling. Results
  are unchanged.
* The colour palette now names every colour used in the figures; the rendered
  colours are unchanged.
* `simulate_biomarkers(n)` and `bootstrap_load(B)` keep their previous
  validation messages. Two `overwrite` arguments now report "must be TRUE or
  FALSE" and `score_imputed()` reports "levels must name 2 distinct existing
  groups".
* The example panel (sbp/hdl/crp) is defined once, in `inst/examples/panel.R`.
* The scripts under `examples/` are no longer part of the distributed tarball;
  they remain available in the repository clone.

## New scoring methods

* `define_protocol(method = "z_mean")` averages risk-oriented z-scores against
  fixed reference means and SDs, with optional `z_cap`; `define_biomarker(transform = "log")`
  scores skewed biomarkers on the log scale.
* `balance = "system"` gives every physiological system equal weight for count
  and z-score methods.
* `bilateral = "split"` divides the tail probability of bilateral quantile
  rules between both tails.

## Normative references

* Protocol `strata` (for example sex) fit thresholds, means and SDs within each
  stratum; unknown strata are rejected. Supplied cutoffs can be stratum-specific
  (`cutoff = c(Female = 50, Male = 40)`).
* Protocol `age` fits growth-chart style references: quantile regression on a
  natural spline of age for count methods (quantreg) and Altman's
  absolute-residual mean/SD model for z-scores. Ages outside the reference
  range receive the status `age_outside_reference` unless `age_extrapolate = TRUE`.
* `plot_centiles()` draws age-specific reference curves with measured values;
  `reference_cutoffs(reference, at = ...)` evaluates them at chosen ages.

## Explaining scores

* `decompose_gap()` splits the difference in mean load between two groups
  exactly into biomarker or system contributions, with bootstrap intervals and a
  waterfall plot.
* `risk_patterns()` tabulates combinations of at-risk biomarkers, their
  co-occurrence (phi, observed vs expected) and the effective number of
  independent indicators; `plot_risk_map()` draws an oncoprint-style map of
  individual contributions.

## Uncertainty and robustness

* `measurement_uncertainty()` propagates user-supplied analytical and
  within-person biological variation (flat or empirical prior) to the exact
  distribution of each participant's score.
* `visit_changes(variation = ...)` classifies risk changes between visits with
  lognormal reference change values.
* `multiverse_grid()` and `multiverse_analysis()` score every common
  operationalization and report participant classification stability,
  agreement and an optional specification curve.
* `bootstrap_load()` also supports z-score protocols and resamples within strata.

## Missing data

* `impute_biomarkers()` performs explicit, auditable imputation: random forests
  within chained equations (mice + randomForest), a native missForest
  implementation, predictive mean matching, CART and Bayesian regression.
* `score_imputed()` pools scores with Rubin's rules (Barnard-Rubin df);
  `fit_reference()` accepts imputations; `check_imputation()` compares methods
  on masked values, including recovery of risk flags.

## Surveys, validation and study design

* Protocol `weights` fit survey-weighted references (weighted quantiles,
  means, SDs and age models); `bootstrap_load(cluster = , design_strata = )`
  resamples primary sampling units within design strata; `survey_summary()`
  reports weighted mean scores and prevalences with cluster-bootstrap intervals.
* `validate_scores()` compares scoring algorithms against binary, continuous or
  survival outcomes with repeated cross-validation (AUC, Brier, Harrell's C,
  R-squared, calibration slope, increment over covariates), refitting
  references inside training folds.
* `medication_rule = "add"` with `medication_effects` applies Tobin-style
  treatment corrections to measured values (count and z-score methods).
* `protocol_template()` lists and rebuilds published panels
  with citations (Seeman et al. 1997).
* `model_trajectory()` fits linear mixed models to score trajectories using
  every visit.
* `impute_biomarkers(delta = )` and `tipping_point()` provide delta-adjusted
  sensitivity analyses for data missing not at random.
* `harmonize_cohorts()` equates biomarker distributions to a target cohort for
  multi-cohort analyses, with before/after diagnostics.
* `detect_outliers()` flags univariate (Tukey or robust z) and multivariate
  (robust Mahalanobis) outliers, optionally winsorizes explicitly and reports
  the impact on scores.

## Graphics and infrastructure

* `score_trajectory(group = ...)` summarizes and plots cohort means by group.
* The sensitivity heatmap orders participants by baseline score and shows
  unavailable pairs in dark gray; a reclassification matrix compares count protocols.
* `summary()` returns a printable `summary.al_result`; references have a
  `print()` method; the recorded package version comes from DESCRIPTION.
* Protocol passports use schema 2; schema 1 passports are upgraded on reading.
* `simulate_biomarkers(covariates = TRUE)` adds age and sex.
* New optional dependencies: quantreg, patchwork, mice and randomForest.

## Core workflow

* Added plot_data() for user-selected ggplot2 geometries, axes and themes.
* Scoring now defaults to graphics = "none"; preset collections are opt-in.
* Simplified the optional theme to a classic ggplot2 style.
* Added automatic boxplots with reproducible individual points and sample sizes,
  plus measured biomarker distributions with fixed reference threshold lines.
* Standardized all package-generated graphics and reports to English;
  the locale argument now accepts only "en". User-supplied labels are preserved.

* Added sensitivity, fixed-reference panel influence and matched score comparisons.
* Added calibration bootstrap with failed-resample reporting.
* Added repeated-visit trajectories and endpoint-aware paired changes.
* Added versioned protocol passports, reference comparability and source metadata.
* Added data-quality audit, explicit unit conversion and threshold graphics.
* Added auditable per-biomarker medication rules.
* Added self-contained HTML reports with aggregate results and software versions.


* Explicit biomarker specifications and reusable reference thresholds.
* Quantile and supplied-cutoff counts with auditable components and coverage.
* Automatic editable ggplot2 graphics and explicit PDF/PNG export.
