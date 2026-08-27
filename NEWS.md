# TaxaTimeImpute 0.99.0

First submission to Bioconductor.

## Features

* Data is supplied in the form a study produces: an abundance table with taxa
  in rows and one column per sample, plus a metadata table saying which
  sample came from which subject at which time. The three metadata columns
  are named through `sample_col`, `subject_col` and `time_col`, so they can
  be called anything. *Subject* refers to whatever the repeated measurements
  were taken on, with no assumption of species.
  `tti_from_metadata()` performs the conversion and reports the design, both
  on screen and in its result. `tti_run()` accepts the same arguments and
  returns the completed table under the caller's own sample names, ordered
  by subject and then time. A column created for a subject-timepoint that
  had no sample is named from its subject and time and flagged in the
  returned metadata, so imputed samples stay distinguishable from measured
  ones.
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
* `cluster_method` chooses how subjects are grouped once a taxon has been
  fitted. `"fpca"`, the default, runs k-means on the FPC scores.
  `"kmeans_fd"` runs functional k-means on the fitted curves through
  **fda.usc**, comparing whole trajectories rather than truncated score
  vectors, and falls back to `"fpca"` on the small sets of curves where
  `kmeans.fd()` cannot fit.
* Input is validated before anything is computed, and each problem reports
  its own cause: sample columns that are not numeric, duplicated taxon
  names, a mask naming subjects or times absent from the table, a mask that
  covers every sample, a design with only one subject, and malformed
  `mask_list` entries. A sample column holding values for some taxa and
  `NA` for others is rejected rather than silently ignored, since the method
  imputes whole missing samples.
* Fits that fall back or fail are reported once at the end of a run, rather
  than per taxon, covering FPCA fits that could not be produced, fits resting
  on three or fewer curves, and functional k-means fallbacks. `tti_metrics()`
  warns when a fit has no observed values to score against.
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
* Imputing a `SummarizedExperiment` that has no row names no longer fails.
  Generated feature names were written onto the completed matrix while the
  object itself still had none, which the assay setter rejected.
