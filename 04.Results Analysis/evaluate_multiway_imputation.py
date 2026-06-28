"""Evaluate imputation results from BAMITA (the ``MultiwayImputation`` tool)
against the ground truth -- the BAMITA counterpart of ``evaluate_deepmicrogen``.

BAMITA writes CLR profiles in exactly the same shape as DeepMicroGen (``#OTU ID``
index, taxa rows x ``subject.timepoint`` columns), so the whole metric machinery
is reused unchanged: this module imports the helpers, the ``COLUMNS`` schema and
the path anchors from ``evaluate_deepmicrogen`` and only redefines how results
are discovered, because BAMITA's output tree has one extra level:

    DeepMicroGen : 04.Benchmarking/<dataset>/DeepMicroGen/<MECH>/<ratio-or-params>/*.csv
    BAMITA       : 04.Benchmarking/<dataset>/MultiwayImputation/<MECH>/<level>/<ratio>/*.csv

BAMITA is not tuned and has no learning-rate/dropout/epoch hyperparameters, so
those three DeepMicroGen-only columns are dropped from the BAMITA schema.

Output: results/MultiwayImputation/<dataset>_MultiwayImputation_<MECH>_performance.csv
"""
import glob
import os

import pandas as pd

from evaluate_deepmicrogen import (
    BENCH_DIR, DATASETS_DIR, RESULTS_DIR, COLUMNS,
    calc_performance, find_mask_file,
)

SEL_TOOL = "MultiwayImputation"        # the R package implementing BAMITA

# DeepMicroGen-only hyperparameter columns -- not applicable to BAMITA
NN_PARAM_COLS = ["Learning Rate", "Dropout Rate", "Number of Epochs"]
COLUMNS_BAMITA = [c for c in COLUMNS if c not in NN_PARAM_COLS]

# BAMITA is run without tuning on every dataset; only ``has_group`` differs.
DATASETS = [
    {"name": "Helminth infection and colon cancer", "has_group": True},
    {"name": "DIABIMMUNE three country",            "has_group": True},
    {"name": "BONUS",                               "has_group": False},
]


def run_dataset(cfg):
    """Evaluate every BAMITA result for one dataset; return a list of records.

    Mirrors ``evaluate_deepmicrogen.run_dataset`` but walks the deeper
    ``<MECH>/<level>/<ratio>/`` tree (so ``missing_method`` is three directories
    up) and omits the DeepMicroGen-only NN hyperparameter columns (BAMITA is not
    tuned). Each record's fields line up with ``COLUMNS_BAMITA``."""
    sel_dataset = cfg["name"]
    has_group = cfg["has_group"]
    records = []

    # one extra "*" vs DeepMicroGen for the <level> directory
    for res in glob.glob(os.path.join(BENCH_DIR, sel_dataset, SEL_TOOL, "*", "*", "*", "*.csv")):
        filename = os.path.basename(res).split('.csv')[0]
        print(filename)
        if filename.endswith('_scaled'):           # evaluate the raw (unscaled) output only
            continue

        seed_val = filename.split('seed_')[1].split('_')[0]
        rep_num = filename.split('_rep_')[1].split('.')[0]
        after = filename.split('relative_abundance_')[1].split('_')
        tax_level = '_'.join(after[:2])
        if has_group:
            group_val = after[2]
            pseudo_str = after[3]
        else:
            group_val = "all"
            pseudo_str = after[2]
        sel_pseudo_count = float(pseudo_str)

        # tree is <MECH>/<level>/<ratio>/file -> mechanism is three levels up
        missing_method = os.path.basename(
            os.path.dirname(os.path.dirname(os.path.dirname(res))))

        imputed_df = pd.read_csv(res, index_col=0)
        gt_group = f"{group_val}_" if has_group else ""
        ground_truth_df = pd.read_csv(
            os.path.join(DATASETS_DIR, sel_dataset, "03.processed", "04.sorted",
                         f"sorted_clr_transformation_relative_abundance_"
                         f"{tax_level}_{gt_group}{pseudo_str}.csv"),
            index_col=0)

        print(res)
        mask_df, mask_ratio = find_mask_file(tax_level, group_val, rep_num, seed_val,
                                             sel_dataset, missing_method, has_group)
        if mask_df is None:
            continue

        performance = calc_performance(mask_df, ground_truth_df, imputed_df, sel_pseudo_count)
        if performance is None:                     # no scorable masked samples
            print(f"  skipped (no masked samples): {os.path.basename(res)}")
            continue

        record = [tax_level, group_val, missing_method, mask_ratio,
                  rep_num, sel_pseudo_count, seed_val]
        record.extend(performance)
        records.append(record)

    return records


if __name__ == '__main__':
    tool_dir = os.path.join(RESULTS_DIR, SEL_TOOL)      # results grouped per tool
    os.makedirs(tool_dir, exist_ok=True)

    for cfg in DATASETS:
        records = run_dataset(cfg)
        df = pd.DataFrame(records, columns=COLUMNS_BAMITA)

        # one CSV per dataset AND per missingness mechanism (MCAR / MAR / MNAR)
        for mechanism, sub in df.groupby("Mechanisms of missingness"):
            out = os.path.join(tool_dir, f'{cfg["name"]}_{SEL_TOOL}_{mechanism}_performance.csv')
            sub.to_csv(out, index=False)
            print(f"wrote {out} ({len(sub)} rows)")
