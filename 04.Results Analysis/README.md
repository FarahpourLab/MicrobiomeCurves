# 05. Results Analysis

Evaluation stage of the CANDi benchmark: it scores the imputation outputs in
`04.Benchmarking/` against the ground truth and turns the scores into summary
tables and plots. Two competitors are handled — **DeepMicroGen** and **BAMITA**
(the `MultiwayImputation` R package).

## What's here

```
05. Results Analysis/
├── evaluate_deepmicrogen.py          # score DeepMicroGen imputations  -> results/DeepMicroGen/*.csv
├── evaluate_multiway_imputation.py   # score BAMITA imputations         -> results/MultiwayImputation/*.csv
├── report_deepmicrogen.py            # tables + plots for DeepMicroGen  (shared machinery)
├── report_multiway_imputation.py     # tables + plots for BAMITA        (imports report_deepmicrogen)
└── results/
    ├── DeepMicroGen/
    │   ├── <dataset>_DeepMicroGen_<MECH>_performance.csv   # one row per imputation run
    │   ├── column_dictionary.xlsx                          # what every CSV column means
    │   └── plot/<dataset>/ ...
    └── MultiwayImputation/
        ├── <dataset>_MultiwayImputation_<MECH>_performance.csv
        ├── column_dictionary.xlsx
        └── plot/<dataset>/ ...
```

`<dataset>` = Helminth / DIABIMMUNE three country / BONUS · `<MECH>` = MCAR / MAR / MNAR.

## Pipeline

1. **Evaluate** (`evaluate_*.py`) — for every masked timepoint, compare the true
   vs imputed composition and write one CSV per (dataset × mechanism). Metrics:
   value errors in CLR / relative-abundance / RA-non-zero space (MSE, RMSE, MAE,
   R², Pearson, Spearman, CCC, NRMSE), compositional distances (Aitchison,
   Bray-Curtis, Jensen-Shannon), structural-zero detection (precision/recall/F1),
   and alpha/beta diversity on stochastic counts (richness, evenness, Shannon
   effective, rarefied Bray-Curtis, Jaccard). See `column_dictionary.xlsx`.
   The DeepMicroGen and BAMITA evaluators share all metric code; BAMITA drops the
   three DeepMicroGen-only hyperparameter columns (it is not tuned).

2. **Report** (`report_*.py`) — load all mechanisms of a dataset together and
   write, under `results/<tool>/plot/<dataset>/`:
   - `summary_by_ratio.csv` — mean/std/median per (mechanism, parameter-set, ratio);
   - `MCAR|MAR|MNAR/<param-set>/<metric>/` — `box`, `violin`, `degradation` plots
     across mask ratios, plus a per-parameter-set `diversity_scatter.png`;
   - `MCAR|MAR|MNAR/.../condition_comparison/` — condition/country curves on one
     axis (datasets with >1 condition only);
   - `mechanism_comparison/<param-set>/<metric>/` — MCAR vs MAR vs MNAR on one axis.

   Only the tuned dataset (Helminth + DeepMicroGen) keeps a single parameter-set
   for plotting (`lr=0.001, dr=0.3, ep=3000`); the others have a single fixed set.
