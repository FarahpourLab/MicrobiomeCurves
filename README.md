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
the row names** and one column per sample, and a sample sheet saying which
sample came from which subject at which time. You name the three metadata
columns, so they can be called whatever your sheet already calls them.

```r
library(TaxaTimeImpute)

counts    # taxa in row names, one column per sample
#>                    S001 S002 S003 S004 S005
#> Bacteroides         412  380  455  501  366
#> Faecalibacterium    198  221  176  205  240
#> Bifidobacterium      87   64  103   58   77
#> Akkermansia          33   41   28   35   44

meta      # one row per sample
#>   SampleID SubjectID Day
#> 1     S001     SUB01   0
#> 2     S002     SUB01   7
#> 3     S003     SUB01  14
#> 4     S004     SUB02   0
#> 5     S005     SUB02   7

run <- tti_run(
    counts, meta,
    sample_col     = "SampleID",
    subject_col    = "SubjectID",
    time_col       = "Day",
    abundance_type = "raw",   # counts; use "clr" if already transformed
    out_dir        = "results",
    K = 1                     # pool all subjects when fitting each taxon;
                              # K = NULL groups similar trajectories instead
                              # and picks how many groups per taxon
)
#> Transforming raw abundances to centred log-ratios.
#> Design: 120 taxa, 58 samples, 20 subjects, 3 time points.
#>   time points, in order: 0, 7, 14
#>   missing: 2 of 60 subject-timepoints (3.3%).
#>     1 with no sample at all: SUB07 at 14
#>     1 whose sample column is entirely NA: SUB13 at 7
#> Drawing uncertainty for 120 taxa ...
#> Wrote 3 file(s) to results: imputed_abundance.tsv, imputation_log.txt,
#>   uncertainty_by_taxon.pdf
```

Nothing is encoded in a column name. `SubjectID` identifies whatever the
repeated measurements were taken on — a participant, a mouse, a plot, a
bioreactor.

**Time** may be numbers, in which case the real spacing is used, or labels
such as `"Baseline"`, `"Week1"`, `"Week4"`. Labels are placed in order at
equal spacing, taken from the factor's levels if it is a factor and from the
row order otherwise. That is reported, because the spacing changes the fit.

**Abundances** may be raw counts or relative abundances
(`abundance_type = "raw"`, CLR-transformed for you with zeros replaced), or
values you have already transformed (`abundance_type = "clr"`).

**`K`** controls how subjects are grouped before a value is imputed. `K = 1`
pools every subject, which is the fastest and the usual starting point.
`K = NULL` clusters similar trajectories per taxon and chooses the number of
clusters by silhouette width, which helps when subjects fall into distinct
response patterns but takes considerably longer.

### What you get back

`run$completed` carries your own sample names, ordered by subject and then
time. Observed values are unchanged. A column is added only for a
subject-timepoint that had no sample; it is named from its subject and time,
and `run$metadata` marks it `imputed = TRUE`.

With `out_dir`, the run writes:

| file | contents |
| --- | --- |
| `imputed_abundance.tsv` | the completed table |
| `imputation_log.txt` | the design, time points in order, every missing subject-timepoint and why, and every warning raised |
| `uncertainty_by_taxon.pdf` | one page per taxon, each imputed value with its 95% interval |

```
MISSING (2 of 60 subject-timepoints)
  SUB07  at time 14  [absent_sample]  no sample was collected
  SUB13  at time 7   [no_data]  sample S038
```

The uncertainty pages are the thing to look at before treating imputed
values as data: a taxon whose gaps rest on very little information shows wide
intervals, and says so on its own page. Turn them off with `plots = FALSE`.
For raster copies, `plot_format = "png"` (or `"both"`) writes one image per
taxon at `dpi`, which defaults to 300; the PDF is vector and is sharp at any
size, so `dpi` does not apply to it.

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
