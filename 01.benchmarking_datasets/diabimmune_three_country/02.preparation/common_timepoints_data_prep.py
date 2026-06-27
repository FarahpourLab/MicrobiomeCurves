"""Build the shared per-country input (metadata + relative abundance) for DIABIMMUNE.

For each country this restricts the data to that country's common timepoints and
to the subjects that actually have them, producing the single shared dataset all
mechanisms (MCAR/MAR/MNAR) consume.

* common timepoints -- the N most frequent ``collection_month`` values within the
  country (same definition as find_common_timepoints.py);
* kept subjects -- those whose real samples cover ALL common months within
  ``+/- TOLERANCE`` months; and, in addition, subjects missing any MAR-driver
  covariate are dropped (via ``analyze_metadata.analyse_df``) so the shared data
  is complete for the amputation;
* abundance -- taken from rel-species-table.csv (positional ``.0-.4``), original
  subjectIDs kept.

Output (per country, into 03.processed/):
    02.relative_abundance/relative_abundance_l7_species_<country>.csv  (#OTU ID + subjectID.tp)
    01.metadata/metadata_l7_species_<country>.csv  (original metadata format,
        all columns, common-timepoint rows)
"""
import pandas as pd
import numpy as np
import os
import importlib.util
from collections import Counter

# this script lives in <dataset>/02.preparation/, so the dataset root is one level up
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL_TABLE = os.path.join(ROOT, "01.input_data", "rel-species-table.csv")
METADATA = os.path.join(ROOT, "01.input_data", "metadata.csv")
METADATA_OUT_DIR = os.path.join(ROOT, "03.processed", "01.metadata")
REL_OUT_DIR = os.path.join(ROOT, "03.processed", "02.relative_abundance")

TIME_COL = "collection_month"
SEPARATION_COL = "country"
TAX_LEVEL = "l7_species"
N_COMMON = 5      # timepoints the abundance table keeps per subject
TOLERANCE = 0     # exact match: a target month counts as present only if the subject has a sample at exactly that collection_month

# reuse the metadata analysis to know which covariates MAR will use, so patients
# missing any of those driver values can be excluded here (one shared, complete dataset)
_spec = importlib.util.spec_from_file_location("am", os.path.join(os.path.dirname(os.path.abspath(__file__)), "analyze_metadata.py"))
am = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(am)


def split_subject_tp(columns):
    parts = columns.to_series().str.rsplit(".", n=1, expand=True)
    parts.columns = ["subject", "tp"]
    parts["tp"] = parts["tp"].astype(int)
    parts["orig"] = list(columns)
    return parts.reset_index(drop=True)


def common_months_for(months_per_subject, n_common):
    freq = Counter()
    for months in months_per_subject:
        freq.update(months)
    freq = pd.Series(freq)
    return sorted(freq.sort_values(ascending=False).head(n_common).index.tolist())


def covers_all(subject_months, common, tol):
    return all(any(abs(mo - t) <= tol for mo in subject_months) for t in common)


def main():
    rel = pd.read_csv(REL_TABLE, index_col=0)
    rel.index.name = "#OTU ID"

    meta = pd.read_csv(METADATA)
    rel_subjects = sorted(set(pd.Series(rel.columns).str.rsplit(".", n=1, expand=True)[0]))
    meta = meta[meta["subjectID"].isin(rel_subjects)]

    subject_country = meta.groupby("subjectID")[SEPARATION_COL].first()
    months_per_subject = (
        meta.dropna(subset=[TIME_COL])
        .groupby("subjectID")[TIME_COL]
        .apply(lambda s: sorted(set(s.astype(int))))
    )

    col_map = split_subject_tp(rel.columns)
    col_map["country"] = col_map["subject"].map(subject_country)

    os.makedirs(METADATA_OUT_DIR, exist_ok=True)
    os.makedirs(REL_OUT_DIR, exist_ok=True)

    for country in sorted(subject_country.dropna().unique()):
        country_subjects = subject_country.index[subject_country == country]
        mps_c = months_per_subject[months_per_subject.index.isin(country_subjects)]

        common = common_months_for(mps_c, N_COMMON)
        qualifying = [s for s in mps_c.index if covers_all(mps_c[s], common, TOLERANCE)]

        sel = (
            col_map[(col_map["country"] == country) & (col_map["subject"].isin(qualifying))]
            .sort_values(["subject", "tp"])
        )

        sel_samples = sel["orig"].tolist()
        rel_country = rel[sel_samples]

        # metadata in the ORIGINAL style (all columns, one row per sample),
        # restricted to the qualifying subjects' samples at the common timepoints
        meta_country = meta[
            meta["subjectID"].isin(qualifying) & meta[TIME_COL].isin(common)
        ].sort_values(["subjectID", TIME_COL])

        # exclude any patient with a missing value in the MAR driver covariates,
        # so the shared dataset is complete for every mechanism (incl. amputation)
        report = am.analyse_df(meta_country)
        mar_covs = report[report["mar_meaningful"] & ~report["time_correlated"]]["feature"].tolist()
        missing_subjects = (
            meta_country.groupby("subjectID")[mar_covs]
            .apply(lambda d: bool(d.isna().any().any()))
        ) if mar_covs else pd.Series(dtype=bool)
        excluded = sorted(missing_subjects.index[missing_subjects]) if len(missing_subjects) else []
        if excluded:
            meta_country = meta_country[~meta_country["subjectID"].isin(excluded)]
            rel_country = rel_country[[c for c in rel_country.columns
                                       if c.rsplit(".", 1)[0] not in excluded]]

        rel_country.to_csv(
            os.path.join(REL_OUT_DIR, f"relative_abundance_{TAX_LEVEL}_{country}.csv"), index=True)
        meta_country.to_csv(
            os.path.join(METADATA_OUT_DIR, f"metadata_{TAX_LEVEL}_{country}.csv"), index=False)

        print(f"[{country}] common months {common} (+/-{TOLERANCE}) -> "
              f"{meta_country['subjectID'].nunique()} subjects, {rel_country.shape[1]} abundance samples"
              f" | excluded (missing MAR covariate): {excluded if excluded else 'none'}")


if __name__ == "__main__":
    main()
