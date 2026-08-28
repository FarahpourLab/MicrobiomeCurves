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

The input is the two tables a study already has. Here they are built from
scratch, so the whole block runs as it stands:

```r
library(TaxaTimeImpute)

taxa <- c(
    "Bacteroides", "Faecalibacterium", "Bifidobacterium", "Akkermansia"
)
subjects <- paste0("SUB", sprintf("%02d", 1:6))
days <- c(0, 7, 14)

# meta: one row per sample, saying which subject it came from and when.
# Any column names will do; you name them in the call below.
meta <- expand.grid(
    SubjectID = subjects, Day = days,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
meta <- meta[order(meta$SubjectID, meta$Day), ]
meta$SampleID <- sprintf("S%03d", seq_len(nrow(meta)))
meta <- meta[, c("SampleID", "SubjectID", "Day")]
rownames(meta) <- NULL

# SUB03 missed the day-14 visit, so that sample was never collected.
meta <- meta[!(meta$SubjectID == "SUB03" & meta$Day == 14), ]

# counts: a taxa-by-samples matrix. Taxa in the ROW NAMES, one column per
# sample, named however the sequencing run named it. Raw counts here.
set.seed(1)
counts <- matrix(
    rpois(length(taxa) * nrow(meta), lambda = c(420, 200, 80, 35)),
    nrow = length(taxa),
    dimnames = list(taxa, meta$SampleID)
)
```

That gives a 4 x 17 abundance matrix and 17 rows of metadata:

```r
counts[, 1:5]
#>                  S001 S002 S003 S004 S005
#> Bacteroides       407  388  413  411  432
#> Faecalibacterium  218  206  188  199  212
#> Bifidobacterium    91   86   74   88   86
#> Akkermansia        37   38   33   39   35

head(meta)
#>   SampleID SubjectID Day
#> 1     S001     SUB01   0
#> 2     S002     SUB01   7
#> 3     S003     SUB01  14
#> 4     S004     SUB02   0
#> 5     S005     SUB02   7
#> 6     S006     SUB02  14
```

Note there are 17 samples for 6 subjects x 3 time points = 18: SUB03 has no
day-14 row at all. That is the gap the package fills.

```r
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
#> Design: 4 taxa, 17 samples, 6 subjects, 3 time points.
#>   time points, in order: 0, 7, 14
#>   missing: 1 of 18 subject-timepoints (5.6%).
#>     1 with no sample at all: SUB03 at 14
#> Detected 1 missing sample(s) across 6 subject(s) and 3 time point(s).
#> Creating 1 empty column(s) for absent samples.
#> Fitting FPCA model over 4 taxa ...
#> Done: 4 of 4 cell(s) imputed.
#> Drawing uncertainty for 4 taxa ...
#> Wrote 3 file(s) to results: imputed_abundance.tsv, imputation_log.txt,
#>   uncertainty_by_taxon.pdf
```

Nothing is encoded in a column name. `SubjectID` identifies whatever the
repeated measurements were taken on — a participant, a mouse, a plot, a
bioreactor.

**Time** may be numbers or labels such as `"Baseline"`, `"Week1"`,
`"Week4"`. Labels take their order from the factor's levels if the column is
a factor and from the row order otherwise, which is reported so you can check
it.

Note that **the model uses the order of the time points, not their values**:
days 0, 7 and 14 are fitted exactly as 0, 1 and 2 would be, so uneven
intervals are not represented in the fitted curve. The values are kept for
reporting, for the axis of the uncertainty pages, and for naming any column
that has to be created.

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
time, with the values on the CLR scale the model works in:

```r
run$completed[, 1:5]
#>              taxon   S001   S002   S003   S004
#> 1      Bacteroides  1.130  1.116  1.258  1.155
#> 2 Faecalibacterium  0.506  0.483  0.471  0.430
#> 3  Bifidobacterium -0.368 -0.391 -0.461 -0.386
#> 4      Akkermansia -1.268 -1.208 -1.269 -1.200
```

Observed values are unchanged. A column is added only for a
subject-timepoint that had no sample; it is named from its subject and time,
and `run$metadata` marks it so it stays distinguishable from a measured one:

```r
run$metadata[run$metadata$imputed, ]
#>      sample subject time time_label imputed
#> 18 SUB03_14   SUB03   14         14    TRUE
```

With `out_dir`, the run writes:

| file | contents |
| --- | --- |
| `imputed_abundance.tsv` | the completed table |
| `imputation_log.txt` | the design, time points in order, every missing subject-timepoint and why, and every warning raised |
| `uncertainty_by_taxon.pdf` | one page per taxon, each imputed value with its 95% interval |

```
MISSING (1 of 18 subject-timepoints)
  SUB03  at time 14  [absent_sample]  no sample was collected
```

Each uncertainty page shows what an imputed value actually rests on: every
subject's trajectory, the fitted curve through them, and the 95% band around
the prediction for the subject whose sample was missing.

With outlier screening on (`use_outliers = TRUE`, the default) the page draws
two fits — one over all subjects, one with the flagged trajectories removed —
because the gap between them is what screening did to that value. With
screening off there is one fit and no outlier vocabulary on the page at all.

One page per taxon, with that taxon's imputed values side by side as facets.
A taxon with more than six spills onto further pages rather than losing any.
When a run would draw a great many pages, it says so before starting.

`tti_uncertainty(run, taxon)` returns the same intervals as a table. Turn the
drawing off with `plots = FALSE`. For raster copies, `plot_format = "png"`
(or `"both"`) writes one image per page at `dpi`, which defaults to 300; the
PDF is vector and is sharp at any size, so `dpi` does not apply to it.

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

The workflow vignette is available in both formats:

```r
vignette("TaxaTimeImpute-workflow", package = "TaxaTimeImpute")      # HTML
vignette("TaxaTimeImpute-workflow-pdf", package = "TaxaTimeImpute")  # PDF
```

Both are generated from the same source, so they cannot drift apart.

If `vignette()` reports that it is not found, the package was installed
without them. Vignettes are not built by default:

```r
remotes::install_github("FarahpourLab/TaxaTimeImpute", build_vignettes = TRUE)
```

## Citation

```r
citation("TaxaTimeImpute")
```

## License

GPL-3.
