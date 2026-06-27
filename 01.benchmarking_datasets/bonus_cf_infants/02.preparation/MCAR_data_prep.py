"""MCAR (Missing Completely At Random) mask generation for BONUS (CF infants).

Single cohort, 7 common timepoints (tp 0..6). The subject x timepoint presence
matrix is rebuilt from the abundance columns and a fraction of cells is masked
uniformly at random (independent of value or time) -> MCAR.

Input : 03.processed/02.relative_abundance/relative_abundance_l7_species.csv
Output: 03.processed/MCAR/mask/l7_species/<ratio>/
        masked_metadata_l7_species.rep_<rep>.seed_<seed>.csv
"""
import pandas as pd
import numpy as np
import os
import glob
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # dataset root (script in 02.preparation)
PROCESSED = os.path.join(ROOT, "03.processed")
INPUT_DIR = os.path.join(PROCESSED, "02.relative_abundance")

MECHANISM = "MCAR"
TAX_LEVEL = "l7_species"

N_REPS = 10
RATIOS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def split_subject_tp(columns):
    """Sample columns look like 'B005.0' -> subject 'B005', timepoint 0."""
    parts = columns.to_series().str.rsplit(".", n=1, expand=True)
    parts.columns = ["subject", "tp"]
    parts["tp"] = parts["tp"].astype(int)
    return parts.reset_index(drop=True)


def build_presence(columns):
    cm = split_subject_tp(columns)
    cm["present"] = 1
    presence = (
        cm.pivot_table(index="subject", columns="tp", values="present", fill_value=0)
        .reset_index()
    )
    first = presence.columns[0]
    rest = sorted(presence.columns[1:], key=lambda x: int(x))
    return presence[[first] + rest]


def random_mask(presence, subject_col, out_dir, file_name, n_reps=N_REPS, ratios=RATIOS):
    meta = presence[[subject_col]]
    mat = presence.drop(columns=[subject_col])
    n_rows, n_cols = mat.shape
    n_cells = n_rows * n_cols

    for r in ratios:
        k = int(n_cells * r)
        seeds = np.random.SeedSequence().generate_state(n_reps).tolist()
        for rep in range(n_reps):
            tmp = mat.copy()
            rng = np.random.default_rng(seeds[rep])
            flat_idx = rng.choice(n_cells, size=k, replace=False)
            row_idx = flat_idx // n_cols
            col_idx = flat_idx % n_cols
            for rr, cc in zip(row_idx, col_idx):
                tmp.iat[rr, cc] = 0
            masked_df = pd.concat([meta, tmp], axis=1)
            os.makedirs(f"{out_dir}/{r}", exist_ok=True)
            masked_df.to_csv(f"{out_dir}/{r}/{file_name}.rep_{rep}.seed_{seeds[rep]}.csv", index=False)


if __name__ == "__main__":
    mask_root = os.path.join(PROCESSED, MECHANISM, "mask", TAX_LEVEL)
    if os.path.isdir(mask_root):
        shutil.rmtree(mask_root)

    for path in sorted(glob.glob(os.path.join(INPUT_DIR, f"relative_abundance_{TAX_LEVEL}.csv"))):
        rel_abundance = pd.read_csv(path, index_col=0)
        presence = build_presence(rel_abundance.columns)
        random_mask(presence, subject_col="subject", out_dir=mask_root,
                    file_name=f"masked_metadata_{TAX_LEVEL}")
        print(f"[{MECHANISM}] presence {presence.shape[0]} subjects x "
              f"{presence.shape[1] - 1} timepoints, {rel_abundance.shape[1]} samples")
