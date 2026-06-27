"""MNAR (Missing Not At Random) mask generation for BONUS (CF infants).

Single cohort, 7 common timepoints. MNAR here means missingness depends on the
timepoint itself (longitudinal dropout). Implemented DIRECTLY rather than with
pyampute: with few discrete timepoints pyampute ties and falls back to MCAR.
Instead each timepoint ``t`` gets a masking weight that grows with ``t`` (later
visits more likely missing), and present cells are sampled proportional to that
weight up to the target ratio.

Input : 03.processed/02.relative_abundance/relative_abundance_l7_species.csv
Output: 03.processed/MNAR/mask/l7_species/<ratio>/
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

MECHANISM = "MNAR"
TAX_LEVEL = "l7_species"

N_REPS = 10
RATIOS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

LATER_MORE_MISSING = True
STEEPNESS = 2.0


def split_subject_tp(columns):
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


def timepoint_weights(timepoints):
    tps = sorted(timepoints)
    n = len(tps)
    weight = {}
    for rank, tp in enumerate(tps):
        frac = rank / (n - 1) if n > 1 else 0.0
        if not LATER_MORE_MISSING:
            frac = 1.0 - frac
        weight[tp] = float(np.exp(STEEPNESS * frac))
    return weight


def mnar_mask(presence, subject_col, out_dir, file_name, n_reps=N_REPS, ratios=RATIOS):
    mat = presence.drop(columns=[subject_col])
    timepoints = list(mat.columns)
    tp_weight = timepoint_weights(timepoints)

    n_rows, n_cols = mat.shape
    base = mat.to_numpy()
    cells = [(i, j) for i in range(n_rows) for j in range(n_cols) if base[i, j] == 1]
    cell_w = np.array([tp_weight[timepoints[j]] for (_, j) in cells], dtype=float)
    cell_w /= cell_w.sum()
    n_cells = len(cells)

    for r in ratios:
        k = min(int(n_cells * r), n_cells)
        seeds = np.random.SeedSequence().generate_state(n_reps).tolist()
        for rep in range(n_reps):
            rng = np.random.default_rng(seeds[rep])
            tmp = mat.copy()
            chosen = rng.choice(n_cells, size=k, replace=False, p=cell_w)
            for idx in chosen:
                i, j = cells[idx]
                tmp.iat[i, j] = 0
            masked_df = pd.concat([presence[[subject_col]], tmp], axis=1)
            os.makedirs(f"{out_dir}/{r}", exist_ok=True)
            masked_df.to_csv(f"{out_dir}/{r}/{file_name}.rep_{rep}.seed_{seeds[rep]}.csv", index=False)


if __name__ == "__main__":
    mask_root = os.path.join(PROCESSED, MECHANISM, "mask", TAX_LEVEL)
    if os.path.isdir(mask_root):
        shutil.rmtree(mask_root)

    for path in sorted(glob.glob(os.path.join(INPUT_DIR, f"relative_abundance_{TAX_LEVEL}.csv"))):
        rel_abundance = pd.read_csv(path, index_col=0)
        presence = build_presence(rel_abundance.columns)
        mnar_mask(presence, subject_col="subject", out_dir=mask_root,
                  file_name=f"masked_metadata_{TAX_LEVEL}")
        print(f"[{MECHANISM}] presence {presence.shape[0]} subjects x "
              f"{presence.shape[1] - 1} timepoints "
              f"(later_more_missing={LATER_MORE_MISSING}, steepness={STEEPNESS})")
