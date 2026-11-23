# Multiple imputation: impute_biomarkers(), score_imputed(), check_imputation()
# and tipping_point().

test_that("multiple imputation keeps observed values and records imputed cells", {
  skip_if_not_installed("mice")
  d <- imputation_data(); attr(d, "al_units") <- c(a = "u")
  set.seed(2); before <- .Random.seed
  imp <- impute_biomarkers(d, imputation_spec(), method = "pmm", m = 3, maxit = 3, predictors = "sex")
  expect_identical(.Random.seed, before)
  expect_length(imp$datasets, 3)
  expect_identical(imp$mask, is.na(as.matrix(d[c("a", "b", "c")])))
  expect_equal(imp$imputed, c("b", "c"))
  for (x in imp$datasets) {
    expect_false(anyNA(x[c("a", "b", "c")]))
    expect_identical(x$b[!is.na(d$b)], d$b[!is.na(d$b)])
    expect_true(all(x$b[is.na(d$b)] %in% d$b))
    expect_identical(attr(x, "al_units"), c(a = "u"))
  }
  expect_identical(impute_biomarkers(d, imputation_spec(), method = "pmm", m = 3, maxit = 3, predictors = "sex")$datasets, imp$datasets)
  expect_true(all(c("distributions", "convergence") %in% names(imp$plots)))
  expect_silent(lapply(imp$plots, ggplot2::ggplot_build))
  expect_output(print(imp), "Biomarker imputation: pmm")
  norm <- impute_biomarkers(d, imputation_spec(), method = "norm", m = 2, maxit = 2)
  expect_true(all(unlist(lapply(norm$datasets, `[[`, "c")) > 0))
  expect_error(impute_biomarkers(transform(d, z = NA), imputation_spec(), predictors = "z"), "complete")
  expect_error(impute_biomarkers(d, imputation_spec(), predictors = "a"), "Biomarkers are imputed")
  expect_error(impute_biomarkers(d, imputation_spec(), method = "missforest", m = 2), "single imputation")
})

test_that("random-forest imputations recover structured values", {
  skip_if_not_installed("mice"); skip_if_not_installed("randomForest")
  d <- imputation_data()
  rf <- impute_biomarkers(d, imputation_spec(), method = "rf", m = 2, maxit = 3)
  expect_false(anyNA(rf$datasets[[2]]$b))
  mf <- impute_biomarkers(d, imputation_spec(), method = "missforest", m = 1, maxit = 5)
  expect_equal(mf$diagnostics$oob$biomarker, c("b", "c"))
  expect_true(all(is.finite(mf$diagnostics$oob$nrmse)))
  imputed <- is.na(d$b)
  truth <- 2 * d$a[imputed]
  expect_lt(mean(abs(mf$datasets[[1]]$b[imputed] - truth)), 0.7 * mean(abs(stats::median(d$b, na.rm = TRUE) - truth)))
  expect_identical(impute_biomarkers(d, imputation_spec(), method = "missforest", m = 1, maxit = 5)$datasets, mf$datasets)
})

test_that("Rubin's rules and pooled references follow their definitions", {
  r <- allostaticR:::.rubin(c(1, 2, 3), c(0.1, 0.1, 0.1), df_complete = 100)
  total <- 0.1 + (4 / 3) * 1; lambda <- (4 / 3) / total
  df_obs <- (101 / 103) * 100 * (1 - lambda); df <- 1 / (lambda^2 / 2 + 1 / df_obs)
  expect_equal(r$estimate, 2); expect_equal(r$std_error, sqrt(total)); expect_equal(r$df, df)
  expect_equal(r$upper, 2 + stats::qt(0.975, df) * sqrt(total))
  single <- allostaticR:::.rubin(2, 0.1, df_complete = 50)
  expect_equal(single$between, 0); expect_equal(single$df, 50)
  skip_if_not_installed("mice")
  d <- imputation_data()
  imp <- impute_biomarkers(d, imputation_spec(), method = "pmm", m = 3, maxit = 2)
  ref <- fit_reference(imp, imputation_spec(), label = "pooled")
  each <- vapply(imp$datasets, function(x) reference_cutoffs(fit_reference(x, imputation_spec()))$upper, numeric(3))
  expect_equal(reference_cutoffs(ref)$upper, rowMeans(each))
  expect_equal(ref$imputation$m, 3)
  pooled <- score_imputed(imp, ref, id = "id", group = "sex")
  q <- vapply(pooled$results, function(x) mean(x$scores$score), numeric(1))
  u <- vapply(pooled$results, function(x) stats::var(x$scores$score) / nrow(d), numeric(1))
  expect_equal(pooled$pooled$estimate[1], mean(q))
  expect_equal(pooled$pooled$std_error[1], sqrt(mean(u) + (4 / 3) * stats::var(q)))
  expect_equal(pooled$pooled$target, c("all", "F", "M", "M minus F"))
  expect_equal(pooled$scores$n_imputed, rowSums(is.na(d[c("a", "b", "c")])))
  expect_true(all(pooled$components$risk_share >= 0 & pooled$components$risk_share <= 1))
  expect_equal(sum(pooled$components$imputed), sum(imp$mask))
  expect_output(print(pooled), "pooled over 3 imputations")
  expect_silent(ggplot2::ggplot_build(pooled$plots$pooled_means))
})

test_that("masking checks compare methods against a median baseline", {
  skip_if_not_installed("mice")
  d <- imputation_data()
  ref <- fit_reference(d, imputation_spec())
  check <- check_imputation(d, imputation_spec(), ref, methods = c("pmm", "norm"), prop = 0.2, m = 2, maxit = 2)
  expect_equal(unique(as.character(check$accuracy$method)), c("median", "pmm", "norm"))
  expect_equal(nrow(check$accuracy), 9)
  b <- check$accuracy[check$accuracy$biomarker == "b", ]
  expect_lt(b$nrmse[b$method == "norm"], b$nrmse[b$method == "median"])
  expect_true(all(check$accuracy$flag_agreement >= 0 & check$accuracy$flag_agreement <= 1))
  expect_silent(lapply(check$plots, ggplot2::ggplot_build))
  expect_output(print(check), "Imputation check")
})

test_that("the imputation check and the scoring share one risk-flag rule", {
  skip_if_not_installed("mice")
  d <- imputation_data()
  spec <- imputation_spec()
  ref <- fit_reference(d, spec)
  prop <- 0.2
  check <- check_imputation(d, spec, ref, methods = "pmm", prop = prop, m = 2, maxit = 2, seed = 1L)
  prep <- allostaticR:::.prepare_imputation(d, spec, character())
  cells <- allostaticR:::.with_seed(1L, function() lapply(seq_along(prep$names), function(j) {
    observed <- which(!prep$mask[, j])
    k <- max(1L, round(prop * length(observed)))
    observed[sample.int(length(observed), k)]
  }))
  params <- allostaticR:::.reference_parameters(ref, d)
  for (j in seq_along(prep$names)) {
    name <- prep$names[j]; b <- spec$biomarkers[[j]]
    truth <- d[[name]][cells[[j]]]
    masked <- d
    masked[[name]][cells[[j]]] <- NA
    imputed <- stats::median(masked[[name]], na.rm = TRUE)
    p <- lapply(params[[j]][c("lower", "upper")], `[`, cells[[j]])
    expected <- mean(allostaticR:::.risk_flag(rep(imputed, length(truth)), p, b$direction, spec$inclusive) ==
      allostaticR:::.risk_flag(truth, p, b$direction, spec$inclusive), na.rm = TRUE)
    row <- check$accuracy[check$accuracy$method == "median" & check$accuracy$biomarker == name, ]
    expect_equal(row$flag_agreement, expected)
  }
})

test_that("score_imputed requires two distinct existing groups", {
  skip_if_not_installed("mice")
  d <- imputation_data()
  imp <- impute_biomarkers(d, imputation_spec(), method = "pmm", m = 2, maxit = 2)
  ref <- fit_reference(imp, imputation_spec(), label = "pooled")
  expect_error(score_imputed(imp, ref, id = "id", group = "sex", levels = c("F", "Z")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
  expect_error(score_imputed(imp, ref, id = "id", group = "sex", levels = c("F", "F")),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
  expect_error(score_imputed(imp, ref, id = "id", group = "sex", levels = "F"),
    "levels must name 2 distinct existing groups.", fixed = TRUE)
})

test_that("delta adjustment shifts only imputed values and finds tipping points", {
  skip_if_not_installed("mice")
  set.seed(6)
  d <- data.frame(g = rep(c("A", "B"), each = 60), a = c(rnorm(60, 10), rnorm(60, 11)), b = rnorm(120, 10))
  d$a[c(5:15, 70:80)] <- NA
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")), method = "z_mean")
  ref <- fit_reference(d, spec)
  base <- impute_biomarkers(d, spec, method = "pmm", m = 2, maxit = 2, predictors = "g")
  shifted <- impute_biomarkers(d, spec, method = "pmm", m = 2, maxit = 2, predictors = "g", delta = c(a = 5))
  miss <- is.na(d$a)
  expect_equal(shifted$datasets[[1]]$a[miss], base$datasets[[1]]$a[miss] + 5)
  expect_equal(shifted$datasets[[1]]$a[!miss], d$a[!miss])
  tp <- tipping_point(d, spec, ref, "a", deltas = c(-6, -3, 0), group = "g", m = 2, maxit = 2, predictors = "g")
  expect_equal(nrow(tp$results), 3)
  expect_true(is.na(tp$tipping_delta) || tp$tipping_delta < 0)
  expect_silent(ggplot2::ggplot_build(tp$plots$tipping))
  expect_error(impute_biomarkers(d, spec, delta = c(zz = 1)), "delta")
})

test_that("impute_biomarkers rejects completely missing biomarkers", {
  d <- data.frame(a = rep(NA_real_, 20), b = rnorm(20))
  spec <- define_protocol(list(define_biomarker("a", "s1", "high", "u"), define_biomarker("b", "s2", "high", "u")))
  expect_error(impute_biomarkers(d, spec, method = "pmm", m = 2), "no observed values")
})
