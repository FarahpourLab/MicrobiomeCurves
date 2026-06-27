"""MAR (Missing At Random) mask generation for DIABIMMUNE.

Missingness is produced by pyampute's ``MultivariateAmputation``: the timepoint
is amputed as a function of a COMBINATION of observed covariates -> MAR. Using
several continuous covariates (not a single discrete one) gives a weighted sum
with many distinct values, so pyampute does NOT fall back to MCAR the way the
single-timepoint MNAR would.

The per-country covariate combination is taken from ``analyze_metadata.py``
(meaningful AND non-time-correlated AND non-binary AND non-outcome) -- e.g.
gest_time, bf_length, num_abx_treatments -- so it stays in sync with the analysis
toggles. A median-fill guard handles any residual NaN before amputation.

Input  : 03.processed/02.relative_abundance/relative_abundance_l7_species_<c>.csv
         03.processed/01.metadata/metadata_l7_species_<c>.csv  (driver-complete)
Output : 03.processed/MAR/mask/l7_species/<ratio>/
         masked_metadata_l7_species_<c>.rep_<rep>.seed_<seed>.csv
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

MECHANISM = "MAR"
TAX_LEVEL = "l7_species"
# shared input used by all mechanisms; the MAR covariates were made complete
# (imputed) in 01.metadata by mechanism_data_prep.py
INPUT_REL = os.path.join(PROCESSED, "02.relative_abundance")
INPUT_META = os.path.join(PROCESSED, "01.metadata")

TIME_COL = "time"            # the amputed (incomplete) variable; weight 0 -> MAR
SUBJECT_COL = "subjectID"

N_REPS = 10
RATIOS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# reuse the metadata analysis so the covariate combination matches the report
_spec = importlib.util.spec_from_file_location("am", os.path.join(os.path.dirname(os.path.abspath(__file__)), "analyze_metadata.py"))
am = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(am)


def mar_covariates(meta_path):
    """Per-country MAR covariate combination: meaningful & not time-correlated."""
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


def build_design(sel, subject_cov, covariates):
    """Per-sample design matrix for pyampute: index = sample name (subjectID.tp),
    columns = the timepoint plus the (broadcast) subject-level covariates."""
    design = pd.DataFrame(index=sel["orig"].values)
    design[TIME_COL] = sel["tp"].values
    for c in covariates:
        design[c] = sel["subject"].map(subject_cov[c]).values
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

    for rel_path in sorted(glob.glob(os.path.join(INPUT_REL, f"relative_abundance_{TAX_LEVEL}_*.csv"))):
        country = os.path.basename(rel_path).split("_")[-1].replace(".csv", "")
        meta_path = os.path.join(INPUT_META, f"metadata_{TAX_LEVEL}_{country}.csv")

        rel_abundance = pd.read_csv(rel_path, index_col=0)
        meta = pd.read_csv(meta_path)

        covariates = mar_covariates(meta_path)
        subject_cov = meta.groupby(SUBJECT_COL)[covariates].first()

        # safety net: amputation needs complete covariates. mechanism_data_prep.py
        # already imputed them in the shared metadata, but median-fill any residual
        # NaN here too so pyampute never sees a missing value.
        if subject_cov.isna().any().any():
            for c in covariates:
                if subject_cov[c].isna().any():
                    subject_cov[c] = subject_cov[c].fillna(subject_cov[c].median())
            print(f"  [warn] {country}: median-filled residual NaN covariate(s) before amputation")

        sel = split_subject_tp(rel_abundance.columns).sort_values(["subject", "tp"])
        presence = build_presence(sel)
        design = build_design(sel, subject_cov, covariates)

        # MAR pattern: missingness of the timepoint depends on the covariates
        # (equal weights); the timepoint itself carries weight 0.
        patterns = [{
            "incomplete_vars": [TIME_COL],
            "weights": {TIME_COL: 0, **{c: 1 for c in covariates}},
            "mechanism": "MAR",
        }]

        amputed_masks(presence, design, patterns, mask_root,
                      file_name=f"masked_metadata_{TAX_LEVEL}_{country}")

        print(f"[{MECHANISM}/{country}] {presence.shape[0]} subjects x "
              f"{presence.shape[1] - 1} timepoints | covariates={covariates}")
