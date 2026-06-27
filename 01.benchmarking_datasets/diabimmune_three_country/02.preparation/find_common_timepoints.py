"""Find the COMMON timepoints (collection_month values) per country for DIABIMMUNE.

The abundance table (rel-species-table.csv) keeps only 5 timepoints per subject,
labelled positionally as ``.0 .. .4``, but the real sampling times live in
metadata.csv -> ``collection_month``. Subjects do not share one identical set of
months (a strict intersection is empty) and some have extra timepoints, so
"common" is defined, within each country, as the ``N`` months that appear in the
most subjects. The script reports those months and, per subject, which common
months are present / missing. It is exploratory -- the actual filtering happens
in common_timepoints_data_prep.py.
"""
import pandas as pd
import numpy as np
import os
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)  # dataset root
REL_TABLE = os.path.join(ROOT, "01.input_data", "rel-species-table.csv")
METADATA = os.path.join(ROOT, "01.input_data", "metadata.csv")

TIME_COL = "collection_month"
SEPARATION_COL = "country"
N_COMMON = 5   # number of timepoints the abundance table keeps per subject


def months_per_subject(meta_subset):
    """subjectID -> sorted unique collection_month values."""
    return (
        meta_subset.dropna(subset=[TIME_COL])
        .groupby("subjectID")[TIME_COL]
        .apply(lambda s: sorted(set(s.astype(int))))
    )


def common_months(mps, n_common):
    """The n_common most frequent months across the given subjects."""
    freq = Counter()
    for months in mps:
        freq.update(months)
    freq = pd.Series(freq).sort_index()
    common = sorted(freq.sort_values(ascending=False).head(n_common).index.tolist())
    return freq, common


def coverage_table(mps, common):
    """Per-subject presence/absence of the common months."""
    rows = []
    for subj, months in mps.items():
        present = [mo for mo in common if mo in months]
        missing = [mo for mo in common if mo not in months]
        rows.append({
            "subjectID": subj,
            "n_common_present": len(present),
            "has_all_common": len(missing) == 0,
            "present_common_months": present,
            "missing_common_months": missing,
            "all_months": months,
        })
    return pd.DataFrame(rows).sort_values(["has_all_common", "subjectID"])


def analyse_country(country, meta_country):
    mps = months_per_subject(meta_country)
    n_subjects = len(mps)
    freq, common = common_months(mps, N_COMMON)

    print(f"\n=========================== {country} ===========================")
    print(f"subjects: {n_subjects}")
    print("collection_month -> #subjects sampled at that month:")
    print(freq.to_string())
    print(f"\n{N_COMMON} COMMON timepoints (collection_month): {common}")
    print("  coverage of each common month:", {mo: int(freq.get(mo, 0)) for mo in common})

    summary = coverage_table(mps, common)
    summary.insert(0, "country", country)
    n_all = int(summary["has_all_common"].sum())
    print(f"subjects with ALL {N_COMMON} common timepoints: {n_all} / {n_subjects}")
    print("distribution of #common-present per subject:")
    print(summary["n_common_present"].value_counts().sort_index().to_string())

    out_csv = os.path.join(BASE, f"common timepoints/common_timepoints_{country}.csv")
    summary.to_csv(out_csv, index=False)
    print(f"wrote -> {out_csv}")
    return summary


def main():
    rel = pd.read_csv(REL_TABLE, index_col=0)
    rel_subjects = sorted(set(pd.Series(rel.columns).str.rsplit(".", n=1, expand=True)[0]))

    meta = pd.read_csv(METADATA)
    meta = meta[meta["subjectID"].isin(rel_subjects)]

    print(meta["subjectID"].unique().shape[0], "subjects in metadata")
    print(len(rel_subjects), "subjects in abundance table")
    print('===========================================================')

    subject_country = meta.groupby("subjectID")[SEPARATION_COL].first()

    all_summaries = []
    for country in sorted(subject_country.dropna().unique()):
        country_subjects = subject_country.index[subject_country == country]
        meta_country = meta[meta["subjectID"].isin(country_subjects)]
        all_summaries.append(analyse_country(country, meta_country))

    # one combined table across countries as well
    combined = pd.concat(all_summaries, ignore_index=True)
    combined_csv = os.path.join(BASE, "common timepoints/common_timepoints_by_country.csv")
    combined.to_csv(combined_csv, index=False)
    print(f"\nwrote combined -> {combined_csv}")


if __name__ == "__main__":
    main()
