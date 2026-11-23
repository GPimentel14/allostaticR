# Graphics: plot_load(), plot_data(), plot_cutoffs(), plot_centiles(),
# plot_risk_map() and theme_allostatic().

test_that("automatic graphics are editable and use observed denominators", {
  d <- data.frame(a = c(10, 0, NA), b = c(5, 6, 4), g = c("A", "A", "B"))
  r <- allostatic_load(d, cutoff_spec(), graphics = "auto", group = "g")
  expect_named(r$plots, c("distribution", "biomarkers", "coverage", "groups", "boxplot", "cutoffs"))
  expect_true(all(vapply(r$plots, inherits, logical(1), "ggplot")))
  expect_equal(r$plots$biomarkers$data$proportion, c(0.5, 2 / 3))
  expect_equal(r$plots$coverage$data$observed, c(2 / 3, 1))
  expect_equal(sum(r$plots$distribution$data$n), 2)
  expect_equal(sum(r$plots$groups$data$proportion), 1)
  expect_silent(lapply(r$plots, ggplot2::ggplot_build))
  expect_null(allostatic_load(d, cutoff_spec(), graphics = "none")$plots)
})

test_that("empty scores preserve coverage and report omissions", {
  r <- allostatic_load(data.frame(a = NA_real_, b = NA_real_), cutoff_spec(), graphics = "auto")
  expect_named(r$plots, "coverage")
  expect_match(attr(r$plots, "omitted")$distribution, "No valid")
  expect_error(plot_load(r), "unavailable")
  expect_silent(ggplot2::ggplot_build(r$plots$coverage))
})

test_that("proportion graphics retain boundary scores and missing groups", {
  d <- data.frame(a = c(0, 10, NA), b = c(6, 5, 4), g = c(NA, "Missing group", "B"))
  r <- allostatic_load(d, cutoff_spec("proportion", 0.5), graphics = "auto", group = "g", locale = "en")
  built <- NULL
  expect_warning(built <- ggplot2::ggplot_build(r$plots$distribution), NA)
  expect_equal(sum(built$data[[1]]$count), 3)
  expect_warning(ggplot2::ggplot_build(r$plots$groups), NA)
  expect_equal(length(unique(r$plots$groups$data$group)), 3)
})

test_that("boxplots retain scored observations, boundaries and caller RNG", {
  d <- data.frame(a = c(0, 10, 10, NA), b = c(6, 5, 6, 4),
    g = c("Control", "Exposed", "Exposed", "Missing panel"))
  set.seed(19); state <- .Random.seed
  r <- allostatic_load(d, cutoff_spec(), graphics = "auto", group = "g")
  built <- ggplot2::ggplot_build(r$plots$boxplot)
  expect_identical(.Random.seed, state)
  expect_equal(sort(built$data[[1]]$x), c(0, 1, 2))
  expect_equal(sort(built$data[[2]]$xmiddle), c(0, 1.5))
  expect_equal(nrow(built$data[[1]]), sum(!is.na(r$scores$score)))
  expect_identical(built$data[[1]]$y,
    ggplot2::ggplot_build(r$plots$boxplot)$data[[1]]$y)
  expect_equal(plot_load(r, "boxplot")$labels$title, "Allostatic load by group")
  expect_equal(plot_load(r, "cutoffs")$labels$title,
    "Biomarker distributions and reference thresholds")
  expect_error(plot_load(r, locale = "es"), "arg")
  ungrouped <- allostatic_load(d, cutoff_spec(), graphics = "auto")
  expect_equal(ungrouped$plots$boxplot$labels$title, "Allostatic load summary")
  partial <- allostatic_load(d, cutoff_spec("proportion", 0.5), graphics = "auto", group = "g")
  expect_equal(sort(ggplot2::ggplot_build(partial$plots$boxplot)$data[[1]]$x),
    sort(partial$scores$score))
})

test_that("users choose graphics and receive unchanged data with aligned metadata", {
  d <- data.frame(a = c(0, 10, NA), b = c(6, 5, 4),
    id = c("p3", "p1", "p2"), g = c("Control", NA, "Exposed"))
  r <- allostatic_load(d, cutoff_spec(), id = "id", group = "g")
  expect_null(r$plots)
  expect_identical(plot_data(r), r$scores)
  expect_identical(plot_data(r, "cutoffs"), r$reference$cutoffs)
  for (type in c("biomarkers", "systems")) {
    x <- plot_data(r, type)
    expect_identical(x$group, d$g[x$row])
    expect_identical(x$id, d$id[x$row])
  }
  expect_true(anyNA(plot_data(r, "biomarkers")$value))
  p <- ggplot2::ggplot(plot_data(r), ggplot2::aes(x = group, y = score)) +
    ggplot2::geom_boxplot(na.rm = TRUE) +
    ggplot2::labs(x = "Study group", y = "Risk score") + ggplot2::theme_bw()
  expect_equal(p$labels$x, "Study group")
  expect_equal(p$labels$y, "Risk score")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("risk maps support counts, z-scores and declared sampling", {
  skip_if_not_installed("patchwork")
  d <- panel_data(60); d$sbp[1:3] <- NA
  r <- allostatic_load(d, panel_spec(), group = "group")
  expect_silent(patchwork::patchworkGrob(plot_risk_map(r)))
  z <- allostatic_load(d, define_protocol(panel_spec()$biomarkers, method = "z_mean"))
  expect_silent(patchwork::patchworkGrob(plot_risk_map(z)))
  set.seed(5); before <- .Random.seed
  sampled <- plot_risk_map(r, max_participants = 20, seed = 2)
  expect_identical(.Random.seed, before)
  expect_match(sampled$patches$annotation$caption, "random sample of 20 of 60")
})

test_that("the package theme and the exported graphics keep their public labels", {
  expect_s3_class(theme_allostatic(), "theme")
  expect_equal(theme_allostatic()$legend.position, "bottom")
  r <- allostatic_load(data.frame(a = c(0, 10), b = c(5, 6)), cutoff_spec(), graphics = "auto")
  expect_equal(r$plots$distribution$labels$title, "Allostatic load distribution")
  expect_equal(r$plots$coverage$labels$title, "Biomarker completeness")
  expect_true(all(vapply(r$plots, function(p) inherits(p, "ggplot"), logical(1))))
})
