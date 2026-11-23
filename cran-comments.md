## Submission notes

* New submission of allostaticR 0.1.0 (first release on CRAN).

* Local `R CMD check --as-cran` with R 4.6.2 (2026-10-29): 0 ERROR,
  0 WARNING and 1 NOTE (`New submission`), the expected note for a first
  submission. The PDF and HTML manuals build correctly.

* System requirements: none. The package needs only R (>= 4.1.0) and the
  packages listed in `Imports` (`ggplot2`, `splines`, `stats`, `utils`); no
  compiler, no external data and no system library is required.

* The packages listed in `Suggests` (`quantreg`, `patchwork`, `mice`,
  `randomForest`, `survival`, `nlme`, `MASS`, `testthat`, `knitr`,
  `rmarkdown`) are accessed only through `requireNamespace()`, so they are
  never needed to install the package. When one of them is missing, the
  corresponding feature reports a clear message and the test suite skips the
  related cases instead of failing.

* All examples, tests and vignettes run completely offline: no example, test
  or vignette downloads data, and no network access is attempted.

* The package does not write any file outside `tempdir()`, and it does not
  modify `options()`, `par()` or the working directory. Files are written only
  when the user supplies an explicit output path.
