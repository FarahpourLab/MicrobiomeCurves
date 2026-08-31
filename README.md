# MicrobiomeCurves

<!-- badges: start -->
[![check-bioc](https://github.com/FarahpourLab/MicrobiomeCurves/actions/workflows/check-bioc.yaml/badge.svg)](https://github.com/FarahpourLab/MicrobiomeCurves/actions/workflows/check-bioc.yaml)
<!-- badges: end -->

MicrobiomeCurves imputes missing time points in longitudinal microbiome data.

A whole sample may be missing because it was not collected or failed quality control. MicrobiomeCurves estimates each missing value by modelling taxon abundance over time as a smooth trajectory and sharing information across subjects with Functional Principal Component Analysis using the PACE formulation. The covariance surface is estimated from all subjects, which supports sparse studies and samples collected at irregular sets of time points.

The package handles missing whole samples. Each sample column must be either fully observed or entirely `NA`. If a column contains observed values for some taxa and `NA` for others, the package reports an error instead of imputing it.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("FarahpourLab/MicrobiomeCurves",
                        build_vignettes = TRUE)
```

The package has no GitHub-only dependencies and needs no system libraries beyond those required by R. Installation does not download anything else.

## Quick start

The input consists of an abundance table and a metadata table. This example builds both from scratch, so the complete block can be run as written:

```r
library(MicrobiomeCurves)

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

This produces a 4 x 17 abundance matrix and a metadata table with 17 rows:

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

Six subjects measured at three time points would give 18 samples, but this dataset has 17. The sample for SUB03 at day 14 is absent, and the package fills that gap.

```r
run <- mc_run(
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

Column names do not need to encode study information. `SubjectID` can identify any subject on which repeated measurements were taken, such as a human participant, mouse, plot, or bioreactor.

**Time** can be numeric or use labels such as `"Baseline"`, `"Week1"`, and `"Week4"`. For a factor, time points follow the factor levels. Otherwise, they follow row order, and the package reports that order so it can be checked.

**The model uses the order of the time points, not their numeric values.** Days 0, 7, and 14 are fitted in the same way as 0, 1, and 2, so unequal intervals are not represented in the fitted trajectory. The original values remain available for reports, uncertainty plot axes, and names assigned to newly created sample columns.

**Abundance type** tells the package whether transformation is needed. With `abundance_type = "raw"`, raw counts or relative abundances are CLR-transformed after zeros are replaced. Use `abundance_type = "clr"` for values that have already been transformed.

**`K`** determines how subjects are grouped before an imputed value is calculated. `K = 1` pools all subjects, runs fastest, and is the usual starting point. `K = NULL` creates a cluster of similar trajectories for each taxon and selects the number of clusters by silhouette width; this can help when subjects show distinct response patterns, but it takes considerably longer.

See the package vignette for a detailed account of how the method works. It walks through the method figure panel by panel.

### What you get back

`run$completed` retains the original sample names and orders samples first by subject and then by time. Its abundance values are on the CLR scale used by the model:

```r
run$completed[, 1:5]
#>              taxon   S001   S002   S003   S004
#> 1      Bacteroides  1.130  1.116  1.258  1.155
#> 2 Faecalibacterium  0.506  0.483  0.471  0.430
#> 3  Bifidobacterium -0.368 -0.391 -0.461 -0.386
#> 4      Akkermansia -1.268 -1.208 -1.269 -1.200
```

Observed values remain unchanged. A new column is added only when a subject has no sample at a time point. Its name is formed from the subject and time, and `run$metadata` marks it as imputed so it remains distinguishable from a measured sample:

```r
run$metadata[run$metadata$imputed, ]
#>      sample subject time time_label imputed
#> 18 SUB03_14   SUB03   14         14    TRUE
```

When `out_dir` is provided, the run writes:

| file | contents |
| --- | --- |
| `imputed_abundance.tsv` | the completed table |
| `imputation_log.txt` | the design, time points in order, every missing subject-timepoint and why, and every warning raised |
| `uncertainty_by_taxon.pdf` | one page per taxon, each imputed value with its 95% interval |

```
MISSING (1 of 18 subject-timepoints)
  SUB03  at time 14  [absent_sample]  no sample was collected
```

Each uncertainty page shows the evidence behind an imputed value. It includes every subject's trajectory, the fitted trajectory, and the 95% band around the prediction for the subject with the missing sample.

With outlier screening enabled through `use_outliers = TRUE`, which is the default, the page shows two fits. One uses all subjects, while the other excludes flagged trajectories, so their difference shows how outlier screening affected the imputed value. With outlier screening disabled, the page shows one fit and uses no outlier screening terminology.

Each page covers one taxon and displays that taxon's imputed values side by side as facets. If a taxon has more than six imputed values, they continue on additional pages so none are omitted. Before producing a large number of pages, the package reports how many it will draw.

`mc_uncertainty(run, taxon)` returns the same intervals in a table. Set `plots = FALSE` to disable drawing. For raster output, `plot_format = "png"` or `"both"` writes one image per page at `dpi`, which defaults to 300; the PDF uses vector graphics and remains sharp at any size, so `dpi` does not affect it.

All written information is also returned in `run$design` and `run$missing`, so no result is available only on screen. To inspect a design without fitting the model, use `mc_from_metadata()` for conversion and reporting alone.

The bundled example is also available through `mc_demo_data()`.

## With SummarizedExperiment

`mc_run()` is a generic with methods for `SummarizedExperiment` and `TreeSummarizedExperiment`. Subject and time are read from `colData`:

```r
library(SummarizedExperiment)

se  <- mc_as_demo_se()
out <- mc_run(se, subject_col = "subject", time_col = "timepoint", K = 1)

assayNames(out)
#> [1] "clr"     "imputed"

metadata(out)$mc_run$missing
```

The input assay remains unchanged. The completed matrix is added as a new assay, keeping observed and imputed values separate. Subject names may contain dots or spaces, and time points do not need to be integers or evenly spaced.

## Two entry points

| Input | Function | Result |
| --- | --- | --- |
| Data with missing samples | `mc_run()` | Completed table. `true_value` is `NA`, since the values were never observed |
| Data for benchmarking | `mc_prepare()` and `mc_fit()` | Values masked on purpose, scored against the truth with `mc_metrics()` |

`mc_metrics()` is not meaningful for a fit from `mc_run()` because a missing value has no observed value for comparison.

## Uncertainty

`mc_ci()` returns analytic intervals based on the FPCA score covariance or returns bootstrap intervals. `mc_plot()` draws a trajectory and its interval:

```r
data(taxa_demo)

prep <- mc_prepare(taxa_demo, "OTU_ID",
                    mask_list = data.frame(rep = "S01", time = 3))
fit  <- mc_fit(prep, K = 1)

mc_plot(fit, species_name = "Taxon01", rep_id = "S01",
         time_id = 3, ci_method = "analytic")
```

## Documentation

```r
vignette("MicrobiomeCurves-workflow", package = "MicrobiomeCurves")
```

For a reading copy outside R, `Rscript tools/render-docs.R` writes the same document to the package root in HTML and PDF formats.

## Citation

```r
citation("MicrobiomeCurves")
```

## License

GPL-3.
