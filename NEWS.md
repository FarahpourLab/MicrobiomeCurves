# TaxaTimeImpute 0.99.0

First submission to Bioconductor.

## Features

* `tti_run()` imputes missing whole samples in a longitudinal microbiome
  table. Missing samples are detected from the data and do not need to be
  supplied as a mask. Both subject-timepoints whose column is entirely `NA`
  and those with no column are covered.
* `tti_run()` is a generic with methods for `data.frame` and for
  `SummarizedExperiment` and `TreeSummarizedExperiment`. The container method
  reads subject and time from `colData`, adds the completed matrix as a new
  assay, and does not modify the input assay.
* `tti_detect_missing()` reports which samples are missing and why, without
  fitting a model.
* `tti_prepare()`, `tti_fit()`, `tti_impute()` and `tti_metrics()` provide the
  benchmarking workflow, in which observed values are masked and then scored.
* `tti_ci()` returns analytic and bootstrap confidence intervals for an
  imputed value. `tti_plot()` draws a trajectory with its interval.
* Fraiman-Muniz functional-depth outlier screening and silhouette-selected
  trajectory clustering are available as options.
* `taxa_demo` provides a small simulated dataset containing both forms of
  missingness. `tti_as_demo_se()` returns the same data as a
  `TreeSummarizedExperiment`.

## Fixes relative to the unreleased predecessor

* Fitting no longer fails when a subject has fewer than two observations.
  Outlier screening returned one flag per subject that FPCA could use, while
  the caller indexed the result against every subject. This raised
  `'names' attribute [n] must be the same length as the vector [m]`.
* Analytic confidence intervals no longer fail when the target subject has no
  remaining observations. An empty index built with `sapply()` returned a list
  instead of an integer vector, which raised
  `invalid subscript type 'list'`. Such a subject now falls back to the
  population mean curve with the prior score variance.
* `tti_plot()` no longer requires **dplyr** and **patchwork** to be attached
  by the user. Both are now imported.
* Fitting no longer resets the caller's random number stream. Results for a
  given `seed` are unchanged.
