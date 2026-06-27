"""Metadata feature audit and MAR-driver selection for DIABIMMUNE.

Audits each per-country metadata file in 03.processed/01.metadata along three
axes and flags which features are usable MAR drivers:

* missingness -- fraction of values (and of patients) that are NA per feature;
* static-ness -- constant overall, constant within a subject (subject-level), or
  time-varying; plus how dominant the most common value is;
* time-correlation -- |Pearson r| of the feature with the timepoint
  (collection_month); a high value means the feature is essentially a time proxy
  and must NOT drive MAR.

Driver-selection filter: ``< 10%`` of patients missing AND not constant AND
``pct_most_frequent < 0.6``. A feature is ``mar_meaningful`` if it additionally is
not time-correlated, not an outcome, not an identifier, and not binary (the
toggles near the top control each rule).

``analyse_df(df)`` returns the per-feature report and is reused by the prep
scripts; running the module writes an Excel workbook
(summary + mar_candidates + per-country sheets).
"""
import pandas as pd
import numpy as np
import os
import glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # dataset root (script in 02.preparation)
META_DIR = os.path.join(ROOT, "03.processed", "01.metadata")
OUT_XLSX = os.path.join(META_DIR, "metadata_feature_report.xlsx")

SUBJECT_COL = "subjectID"
TIME_REF = "collection_month"        # canonical timepoint variable for time-correlation
TIME_CORR_THRESHOLD = 0.5            # |r| >= this => treated as time-correlated
USER_FILTER_TOP = 0.6               # pct_most_frequent must be below this
MAX_PCT_SUBJECTS_MISSING = 0.10     # keep a feature if < this fraction of patients have any missing value
REQUIRE_TIME_VARYING = False        # if True, drop features that are constant within a subject (subject-static)
EXCLUDE_TIME_CORRELATED = False     # if True, time-correlated features are dropped from mar_meaningful
EXCLUDE_BINARY = True               # if True, binary (2-value) features are dropped from mar_meaningful

# outcome variables that should not drive missingness (would confound analysis)
OUTCOME_COLS = {
    "seroconverted", "num_aabs", "totalige", "totalige_log", "totalige_high",
    "allergy_milk", "allergy_egg", "allergy_peanut", "allergy_dustmite",
    "allergy_cat", "allergy_dog", "allergy_birch", "allergy_timothy",
}


def to_numeric(series):
    if series.dtype == bool:
        return series.astype(float)
    return pd.to_numeric(series, errors="coerce")


def classify(n_nonnull, n_unique, is_constant, constant_within_subject):
    if n_nonnull == 0:
        return "all-missing"
    if n_unique == n_nonnull and n_nonnull > 1:
        return "identifier/unique"
    if is_constant:
        return "constant"
    if constant_within_subject:
        return "subject-static"
    return "time-varying"


def analyse_file(path):
    return analyse_df(pd.read_csv(path))


def analyse_df(df):
    n_rows = len(df)
    by_subject = df.groupby(SUBJECT_COL) if SUBJECT_COL in df.columns else None
    n_subjects = df[SUBJECT_COL].nunique() if SUBJECT_COL in df.columns else n_rows
    time_ref = to_numeric(df[TIME_REF]) if TIME_REF in df.columns else None

    rows = []
    for col in df.columns:
        s = df[col]
        n_nonnull = int(s.notna().sum())
        n_missing = int(s.isna().sum())
        # how many distinct subjects have at least one missing value in this column
        if by_subject is not None:
            n_subjects_missing = int(by_subject[col].apply(lambda x: x.isna().any()).sum())
        else:
            n_subjects_missing = int(n_missing > 0)
        n_unique = int(s.nunique(dropna=True))
        is_constant = n_unique <= 1
        pct_top = round(s.value_counts(dropna=True).iloc[0] / n_nonnull, 3) if n_nonnull else np.nan

        if by_subject is not None and col != SUBJECT_COL:
            varies = by_subject[col].nunique(dropna=True) > 1
            constant_within_subject = not bool(varies.any())
            pct_subjects_varying = round(float(varies.mean()), 3)
        else:
            constant_within_subject = is_constant
            pct_subjects_varying = 0.0

        # |correlation| with the timepoint
        if col == TIME_REF:
            corr_time = 1.0
        elif time_ref is not None:
            s_num = to_numeric(s)
            if s_num.notna().sum() > 1 and s_num.nunique() > 1:
                corr_time = abs(s_num.corr(time_ref))
            else:
                corr_time = np.nan
        else:
            corr_time = np.nan

        rows.append({
            "feature": col,
            "dtype": str(s.dtype),
            "n_rows": n_rows,
            "n_missing": n_missing,
            "pct_missing": round(n_missing / n_rows, 3) if n_rows else np.nan,
            "n_subjects_missing": n_subjects_missing,
            "pct_subjects_missing": round(n_subjects_missing / n_subjects, 3) if n_subjects else np.nan,
            "n_unique": n_unique,
            "is_constant": is_constant,
            "pct_most_frequent": pct_top,
            "constant_within_subject": constant_within_subject,
            "pct_subjects_timevarying": pct_subjects_varying,
            "corr_with_time": round(corr_time, 3) if pd.notna(corr_time) else np.nan,
            "category": classify(n_nonnull, n_unique, is_constant, constant_within_subject),
        })

    report = pd.DataFrame(rows)

    # --- selection flags -----------------------------------------------------
    report["is_outcome"] = report["feature"].isin(OUTCOME_COLS)
    report["is_binary"] = report["n_unique"] == 2
    passes = (
        (report["pct_subjects_missing"] < MAX_PCT_SUBJECTS_MISSING)
        & (~report["is_constant"])
        & (report["pct_most_frequent"] < USER_FILTER_TOP)
    )
    if REQUIRE_TIME_VARYING:
        passes &= ~report["constant_within_subject"]
    report["passes_user_filter"] = passes
    report["time_correlated"] = report["corr_with_time"] >= TIME_CORR_THRESHOLD
    meaningful = (
        report["passes_user_filter"]
        & report["corr_with_time"].notna()                  # numeric & time-assessable
        & (~report["is_outcome"])
        & (report["feature"] != SUBJECT_COL)
        & (~report["category"].isin(["identifier/unique"]))
    )
    if EXCLUDE_TIME_CORRELATED:
        meaningful &= ~report["time_correlated"]
    if EXCLUDE_BINARY:
        meaningful &= ~report["is_binary"]
    report["mar_meaningful"] = meaningful
    return report


def main():
    paths = sorted(glob.glob(os.path.join(META_DIR, "metadata_*.csv")))
    if not paths:
        print(f"No metadata files found in {META_DIR}")
        return

    per_country, summaries = {}, []
    for path in paths:
        country = os.path.basename(path).replace("metadata_", "").replace(".csv", "").split("_")[-1]
        report = analyse_file(path)
        per_country[country] = report
        s = report.copy()
        s.insert(0, "country", country)
        summaries.append(s)
        print(f"[{country}] {report.shape[0]} features | "
              f"pass user filter: {int(report['passes_user_filter'].sum())} | "
              f"meaningful (non-time, non-outcome): {int(report['mar_meaningful'].sum())}")

    combined = pd.concat(summaries, ignore_index=True)

    # features passing the user filter, side by side across countries
    cand = combined[combined["passes_user_filter"]][
        ["country", "feature", "n_subjects_missing", "pct_subjects_missing",
         "pct_most_frequent", "corr_with_time",
         "time_correlated", "is_binary", "is_outcome", "mar_meaningful"]
    ].sort_values(["feature", "country"])

    def write_workbook(path):
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            combined.to_excel(writer, sheet_name="summary", index=False)
            cand.to_excel(writer, sheet_name="mar_candidates", index=False)
            for country, report in per_country.items():
                report.to_excel(writer, sheet_name=country[:31], index=False)

    try:
        write_workbook(OUT_XLSX)
        print(f"\nwrote -> {OUT_XLSX}")
    except PermissionError:
        alt = OUT_XLSX.replace(".xlsx", "_new.xlsx")
        write_workbook(alt)
        print(f"\n{OUT_XLSX} was locked (open in Excel); wrote -> {alt} instead")


if __name__ == "__main__":
    main()
