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

```r
library(TaxaTimeImpute)
data(taxa_demo)

# report what is missing
tti_detect_missing(taxa_demo, taxon_col = "OTU_ID")$missing
#>    rep time        reason    col
#> 1  S02    1        all_na S02.1
#> 2  S04    2 absent_column S04.2
#> 3  S07    4        all_na S07.4

# impute
run <- tti_run(taxa_demo, taxon_col = "OTU_ID", K = 1)
run
#> TaxaTimeImpute run
#>   taxa             : 12
#>   subjects         : 10
#>   missing samples  : 3
#>   cells imputed    : 36 of 36

head(run$completed)
```

Observed values are not modified. Existing columns keep their order and their
values. A column is added only when a subject-timepoint had no sample.

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
