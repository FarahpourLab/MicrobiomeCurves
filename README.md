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

Give it the two tables a study already has: an abundance table with **taxa in
the row names** and one column per sample, and a metadata sheet saying which
sample came from which subject at which time. You name the three metadata
columns, so they can be called anything.

```r
library(TaxaTimeImpute)

counts    # taxa in row names, one column per sample
#>             RUN_0031 RUN_0044 RUN_0052 RUN_0067
#> Bacteroides       12       15       10       13
#> Prevotella         4        6        5        7
#> Akkermansia        9       11        8       10

meta      # one row per sample
#>    library animal day
#> 1 RUN_0031    M01   0
#> 2 RUN_0044    M01   7
#> 3 RUN_0052    M02   0
#> 4 RUN_0067    M02   7

run <- tti_run(
    counts, meta,
    sample_col     = "library",
    subject_col    = "animal",
    time_col       = "day",
    abundance_type = "raw",       # counts; use "clr" if already transformed
    out_dir        = "results",
    K = 1
)
#> Transforming raw abundances to centred log-ratios.
#> Design: 8 taxa, 23 samples, 6 subjects, 4 time points.
#>   time points, in order: baseline (1), week1 (2), week4 (3), week8 (4)
#>   missing: 2 of 24 subject-timepoints (8.3%).
#>     1 with no sample at all: M03 at week4
#>     1 whose sample column is entirely NA: M04 at week4
#> Wrote imputed_abundance.tsv and imputation_log.txt to results
```

Nothing is encoded in a column name. `subject` is whatever the repeated
measurements were taken on — a mouse, a participant, a plot, a bioreactor.

**Time** may be numbers, in which case the real spacing is used, or labels
such as `"baseline"`, `"week1"`, `"week4"`. Labels are placed in order at
equal spacing, taken from the factor's levels if it is a factor and from the
row order otherwise. That is reported, because the spacing changes the fit.

**Abundances** may be raw counts or relative abundances
(`abundance_type = "raw"`, CLR-transformed for you, zeros replaced), or
values you have already transformed (`abundance_type = "clr"`).

### What you get back

`run$completed` carries your own sample names, ordered by subject and then
time. Observed values are unchanged. A column is added only for a
subject-timepoint that had no sample; it is named from its subject and time,
and `run$metadata` marks it `imputed = TRUE`.

With `out_dir`, two files are written:

- `imputed_abundance.tsv` — the completed table
- `imputation_log.txt` — the time points in order, the counts of samples,
  subjects and time points, every subject-timepoint that was missing and
  why, and every warning raised

```
MISSING (2 of 24 subject-timepoints)
  M03  at time week4  [absent_sample]  no sample was collected
  M04  at time week4  [no_data]  sample RUN_0016
```

Everything written is also returned, in `run$design` and `run$missing`, so
nothing is available only on screen. To inspect a design without fitting,
`tti_from_metadata()` does the conversion and reporting alone.

The bundled example is available in this form with `tti_demo_data()`.

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
