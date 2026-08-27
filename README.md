# TaxaTimeImpute

<!-- badges: start -->
[![check-bioc](https://github.com/FarahpourLab/TaxaTimeImpute/actions/workflows/check-bioc.yaml/badge.svg)](https://github.com/FarahpourLab/TaxaTimeImpute/actions/workflows/check-bioc.yaml)
<!-- badges: end -->

Imputation of missing timepoints in longitudinal microbiome data.

Longitudinal microbiome studies often lose whole samples. A subject misses a
visit, or a library fails quality control, and the affected subject-timepoint
is absent from the abundance table. TaxaTimeImpute estimates the missing
values by modelling each taxon's abundance over time as a smooth curve and
sharing information across subjects using Functional Principal Component
Analysis (FPCA, PACE formulation). Because the covariance surface is estimated
from all subjects together, sparse and irregularly sampled designs are
supported.

The package handles missing whole samples. A sample column must therefore be
either fully observed or entirely `NA`; a column holding values for some taxa
and `NA` for others is a different problem, and is reported as an error rather
than imputed.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("FarahpourLab/TaxaTimeImpute")
```

There are no GitHub-only dependencies, no system libraries beyond those R
itself requires, and nothing is downloaded during installation.

## Quick start

Give it the two tables a study already has: an abundance table with taxa in
rows and samples in columns, and a metadata sheet saying which sample came
from which subject at which time. You name the three metadata columns, so
they can be called anything.

```r
library(TaxaTimeImpute)

counts    # taxa in rows, one column per sample
#>         taxon M01_d0 M01_d7 M02_d0 M02_d7
#> 1 Bacteroides   1.20   1.51   0.98   1.14
#> 2  Prevotella   0.41   0.62   0.55   0.47

meta      # one row per sample
#>   sample subject day
#> 1 M01_d0     M01   0
#> 2 M01_d7     M01   7
#> 3 M02_d0     M02   0
#> 4 M02_d7     M02   7

run <- tti_run(
    counts,
    metadata    = meta,
    sample_col  = "sample",
    subject_col = "subject",
    time_col    = "day",
    taxon_col   = "taxon",
    K = 1
)
#> Design: 8 taxa, 23 samples, 6 subjects, 4 time points.
#>   time points: 0, 7, 14, 21
#>   missing: 2 of 24 subject-timepoints (8.3%).
#>     1 with no sample at all: M03 at 14
#>     1 whose sample column is entirely NA: M05 at 21

head(run$completed)
```

`subject` is whatever the repeated measurements were taken on — a mouse, a
participant, a plot, a bioreactor. Nothing in the package assumes a species.

The completed table comes back under your own sample names, ordered by
subject and then time. Observed values are not modified. A column is added
only for a subject-timepoint that had no sample at all; it is named from its
subject and time, and `run$metadata` marks it `imputed = TRUE` so it stays
distinguishable from a measured sample.

Everything printed above is also returned, in `run$design`, so nothing is
available only by reading the console. To inspect a design without fitting
anything, `tti_from_metadata()` does the conversion and reporting on its own.

If your table already encodes subject and time in its column names as
`"<subject>.<time>"`, omit `metadata` and the original interface still
applies:

```r
data(taxa_demo)
run <- tti_run(taxa_demo, taxon_col = "OTU_ID", K = 1)
```

## With SummarizedExperiment

`tti_run()` is a generic with a method for `SummarizedExperiment` and
`TreeSummarizedExperiment`. Subject and time are read from `colData`:

```r
library(SummarizedExperiment)

se  <- tti_as_demo_se()
out <- tti_run(se, subject_col = "subject", time_col = "timepoint", K = 1)

assayNames(out)
#> [1] "clr"     "imputed"

metadata(out)$tti_run$missing
```

The input assay is not modified. The completed matrix is added as a new assay,
which keeps observed and imputed values separate. Subject names may contain
dots or spaces, and time points need not be integers or evenly spaced.

## Two entry points

| Input | Function | Result |
| --- | --- | --- |
| Data with missing samples | `tti_run()` | Completed table. `true_value` is `NA`, since the values were never observed |
| Data for benchmarking | `tti_prepare()` and `tti_fit()` | Values masked on purpose, scored against the truth with `tti_metrics()` |

`tti_metrics()` is not meaningful for a `tti_run()` fit, because the missing
values have no observed counterpart.

## Uncertainty

`tti_ci()` returns analytic intervals from the FPCA score covariance, or
bootstrap intervals. `tti_plot()` draws a trajectory with its interval:

```r
prep <- tti_prepare(taxa_demo, "OTU_ID",
                    mask_list = data.frame(rep = "S01", time = 3))
fit  <- tti_fit(prep, K = 1)

tti_plot(fit, species_name = "Taxon01", rep_id = "S01",
         time_id = 3, ci_method = "analytic")
```

## Documentation

```r
vignette("TaxaTimeImpute-workflow", package = "TaxaTimeImpute")
```

## Citation

```r
citation("TaxaTimeImpute")
```

## License

GPL-3.
