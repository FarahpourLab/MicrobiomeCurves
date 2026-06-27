"""MAR (Missing At Random) mask generation for BONUS (CF infants).

Single cohort, 7 common timepoints. Missingness of the timepoint depends on
observed covariate(s) -> MAR, produced by pyampute. The covariate set comes from
analyze_metadata.py (meaningful drivers, excluding time proxies and outcomes) and
is read per-sample from the metadata (BONUS covariates are time-varying). A
continuous covariate such as ``weight %`` has many distinct values, so pyampute
does not fall back to MCAR.

Input : 03.processed/02.relative_abundance/relative_abundance_l7_species.csv
        03.processed/01.metadata/metadata_l7_species.csv
Output: 03.processed/MAR/mask/l7_species/<ratio>/
        masked_metadata_l7_species.rep_<rep>.seed_<seed>.csv
"""
import pandas as pd
import numpy as np
import os
import glob
import shutil
import importlib.util
from pyampute.ampute import MultivariateAmputation

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # dataset root (script in 02.preparation)
PROCESSED = os.path.join(ROOT, "03.processed")
INPUT_REL = os.path.join(PROCESSED, "02.relative_abundance")
INPUT_META = os.path.join(PROCESSED, "01.metadata")

MECHANISM = "MAR"
TAX_LEVEL = "l7_species"
TIME_COL = "time"            # the amputed (incomplete) variable; weight 0 -> MAR
SAMPLE_COL = "sample"        # metadata column matching the abundance column names

N_REPS = 10
RATIOS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

_spec = importlib.util.spec_from_file_location("am", os.path.join(os.path.dirname(os.path.abspath(__file__)), "analyze_metadata.py"))
am = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(am)


def mar_covariates(meta_path):
    """Meaningful & not time-correlated covariate combination from the analysis."""
    r = am.analyse_file(meta_path)
    return r[r["mar_meaningful"] & ~r["time_correlated"]]["feature"].tolist()


def split_subject_tp(columns):
    parts = columns.to_series().str.rsplit(".", n=1, expand=True)
    parts.columns = ["subject", "tp"]
    parts["tp"] = parts["tp"].astype(int)
    parts["orig"] = list(columns)
    return parts.reset_index(drop=True)


def build_presence(sel):
    sel = sel.copy()
    sel["present"] = 1
    presence = (
        sel.pivot_table(index="subject", columns="tp", values="present", fill_value=0)
        .reset_index()
    )
    first = presence.columns[0]
    rest = sorted(presence.columns[1:], key=lambda x: int(x))
    return presence[[first] + rest]


def build_design(rel_columns, meta, covariates):
    """Per-sample design matrix for pyampute: index = sample (B005.0), columns =
    the timepoint plus the per-sample covariate values from the metadata."""
    cm = split_subject_tp(rel_columns)                 # subject, tp, orig(sample)
    design = cm.set_index("orig")[["tp"]].rename(columns={"tp": TIME_COL})
    cov = meta.set_index(SAMPLE_COL)[covariates]
    design = design.join(cov, how="left")

    # amputation needs complete covariates: median-fill any residual NaN
    for c in covariates:
        if design[c].isna().any():
            design[c] = design[c].fillna(design[c].median())
    return design


def amputed_masks(presence, design, patterns, out_dir, file_name):
    for r in RATIOS:
        seeds = np.random.SeedSequence().generate_state(N_REPS).tolist()
        for rep in range(N_REPS):
            masked_df = presence.copy()
            ma = MultivariateAmputation(prop=r, patterns=patterns, seed=seeds[rep])
            amputed = ma.fit_transform(design)

            na_samples = amputed.loc[amputed[TIME_COL].isna()].index
            for sample_name in na_samples:
                subject, tp = sample_name.rsplit(".", 1)
                row = masked_df.index[masked_df["subject"] == subject]
                masked_df.loc[row[0], int(tp)] = 0

            os.makedirs(f"{out_dir}/{r}", exist_ok=True)
            masked_df.to_csv(f"{out_dir}/{r}/{file_name}.rep_{rep}.seed_{seeds[rep]}.csv", index=False)


if __name__ == "__main__":
    mask_root = os.path.join(PROCESSED, MECHANISM, "mask", TAX_LEVEL)
    if os.path.isdir(mask_root):
        shutil.rmtree(mask_root)

    for rel_path in sorted(glob.glob(os.path.join(INPUT_REL, f"relative_abundance_{TAX_LEVEL}.csv"))):
        meta_path = os.path.join(INPUT_META, f"metadata_{TAX_LEVEL}.csv")
        rel_abundance = pd.read_csv(rel_path, index_col=0)
        meta = pd.read_csv(meta_path)

        covariates = mar_covariates(meta_path)
        if not covariates:
            raise SystemExit("no meaningful, non-time-correlated MAR covariates found")

        sel = split_subject_tp(rel_abundance.columns).sort_values(["subject", "tp"])
        presence = build_presence(sel)
        design = build_design(rel_abundance.columns, meta, covariates)

        patterns = [{
            "incomplete_vars": [TIME_COL],
            "weights": {TIME_COL: 0, **{c: 1 for c in covariates}},
            "mechanism": "MAR",
        }]

        amputed_masks(presence, design, patterns, mask_root,
                      file_name=f"masked_metadata_{TAX_LEVEL}")
        print(f"[{MECHANISM}] presence {presence.shape[0]} subjects x "
              f"{presence.shape[1] - 1} timepoints | covariates={covariates}")
