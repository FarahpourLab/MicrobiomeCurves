# TIM-FDA
Time series Imputation for Microbiome using Functional Data Analysis


Time series Imputation for Microbiome using Functional Data Analysis: a benchmark
for **imputing missing microbiome timepoints** in longitudinal microbiome
profiles. Datasets are masked under three missingness mechanisms (MCAR / MAR /
MNAR), competitor methods impute the masked timepoints, and the imputations are
scored against the held-out truth.

The pipeline runs top-to-bottom through the numbered folders. Each folder has its
own README with the details.

```
.
├── 01.benchmarking_datasets/   # the datasets: raw inputs, preparation, processed data + masks
├── 02.run_competitors/         # scripts that run each competitor on each dataset
├── 03.imputation_results/      # the imputed outputs produced by the competitors
├── 04.Results Analysis/        # scoring the imputations into metrics, tables and plots
```

See each folder's README for what it contains and how to use it.
