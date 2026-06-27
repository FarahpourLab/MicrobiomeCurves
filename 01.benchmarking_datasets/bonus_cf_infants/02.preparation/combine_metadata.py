"""Build the shared input (metadata + relative abundance) for BONUS (CF infants).

Combines the two sheets of the supplementary Excel into one per-sample metadata
table and aligns it with the species profile:

* SuppTable1 = patient-level metadata -> broadcast (constant) across the
  patient's timepoints;
* SuppTable2 = sample-level (per-month) metadata.

Samples are matched to profile_Species.csv and renamed to the profile convention
``<participant>.<tp>`` (tp 0..6). The profile's tp index maps to a FIXED month
grid (the 7 most frequent months: 3, 4, 5, 6, 8, 10, 12). A patient is kept only
if it has ALL 7 common timepoints; subjects missing any are EXCLUDED from both
the metadata and the abundance (off-grid extra months are reported and dropped).

Output (into 03.processed/):
    01.metadata/metadata_l7_species.csv
    02.relative_abundance/relative_abundance_l7_species.csv
"""
import os
import re
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))            # 02.preparation
INPUT_DIR = os.path.join(BASE, "..", "01.input_data")
EXCEL = os.path.join(INPUT_DIR, "41591_2019_714_MOESM1_ESM.xlsx")
PROFILE = os.path.join(INPUT_DIR, "profile_Species.csv")
OUT_META_DIR = os.path.join(BASE, "..", "03.processed", "01.metadata")
OUT_REL_DIR = os.path.join(BASE, "..", "03.processed", "02.relative_abundance")

HEADER_ROW = 4   # the two sheets have title/description rows above the header


def main():
    # --- profile: valid samples = <participant>.<tp> ------------------------
    prof_cols = pd.read_csv(PROFILE, nrows=0).columns.tolist()[1:]
    prof = pd.DataFrame({"sample": prof_cols})
    parts = prof["sample"].str.rsplit(".", n=1, expand=True)
    prof["participant"] = parts[0]
    prof["tp"] = parts[1].astype(int)
    prof_participants = set(prof["participant"])

    # --- patient-level (SuppTable1) -----------------------------------------
    s1 = pd.read_excel(EXCEL, sheet_name="SuppTable1", skiprows=HEADER_ROW)
    s1 = s1[s1["participant_id"].astype(str).str.match(r"^B\d+$", na=False)].copy()

    # --- sample-level (SuppTable2) -------------------------------------------
    s2 = pd.read_excel(EXCEL, sheet_name="SuppTable2", skiprows=HEADER_ROW)
    s2 = s2[s2["sample"].astype(str).str.match(r"^B\d+-M\d+$", na=False)].copy()
    s2["participant"] = s2["sample"].str.extract(r"^([^-]+)-M")
    s2["month"] = s2["sample"].str.extract(r"-M(\d+)$")[0].astype(int)
    # keep the original metadata sample id but avoid colliding with the profile
    # sample name ("<participant>.<tp>")
    s2 = s2.rename(columns={"sample": "orig_sample"})

    # --- fixed month grid (7 most frequent months) --------------------------
    common = s2[s2["participant"].isin(prof_participants)]
    grid = sorted(common["month"].value_counts().head(7).index.tolist())
    if len(grid) != 7:
        raise SystemExit(f"expected a 7-month grid, got {grid}")
    month_to_tp = {m: i for i, m in enumerate(grid)}
    print(f"month grid (tp 0..6): {grid}")

    # --- report patients with more timepoints than the grid -----------------
    print("\nPatients with EXTRA (off-grid) timepoints in metadata:")
    n_flagged = 0
    for pid, g in common.groupby("participant"):
        extra = sorted(set(g["month"]) - set(grid))
        if extra:
            n_flagged += 1
            print(f"  {pid}: {len(g)} metadata months {sorted(g['month'].tolist())} "
                  f"-> extra (dropped): {extra}")
    print(f"({n_flagged} patient(s) flagged)\n")

    # --- map each grid sample to its profile timepoint ----------------------
    s2_grid = s2[s2["month"].isin(grid)].copy()
    s2_grid["tp"] = s2_grid["month"].map(month_to_tp)
    s2_grid = s2_grid.drop_duplicates(["participant", "tp"], keep="first")

    # --- assemble: every profile sample gets sample- + patient-level fields --
    combined = prof.merge(s2_grid, on=["participant", "tp"], how="left")
    combined = combined.merge(s1, left_on="participant", right_on="participant_id", how="left")

    # column order: sample, participant, tp/month, patient-level, sample-level
    patient_cols = [c for c in s1.columns if c != "participant_id"]
    sample_cols = [c for c in s2.columns if c not in ("orig_sample", "participant", "month")]
    ordered = (["sample", "participant", "tp", "month", "orig_sample"]
               + ["participant_id"] + patient_cols + sample_cols)
    combined = combined[[c for c in ordered if c in combined.columns]]
    combined = combined.sort_values(["participant", "tp"]).reset_index(drop=True)

    # Keep only subjects that have ALL common (grid) timepoints. A subject with a
    # padded timepoint (missing grid month -> NaN orig_sample) is excluded
    # entirely, from both the metadata and the abundance.
    complete = combined.groupby("participant")["orig_sample"].apply(lambda s: bool(s.notna().all()))
    kept = set(complete.index[complete])
    excluded = sorted(complete.index[~complete])

    meta_out = combined[combined["participant"].isin(kept)].reset_index(drop=True)

    # filter the abundance profile to the kept subjects' samples
    prof_ab = pd.read_csv(PROFILE, index_col=0)
    prof_ab.index.name = "#OTU ID"
    keep_cols = [c for c in prof_ab.columns if c.rsplit(".", 1)[0] in kept]
    ab_out = prof_ab[keep_cols]

    os.makedirs(OUT_META_DIR, exist_ok=True)
    os.makedirs(OUT_REL_DIR, exist_ok=True)
    meta_path = os.path.join(OUT_META_DIR, "metadata_l7_species.csv")
    rel_path = os.path.join(OUT_REL_DIR, "relative_abundance_l7_species.csv")
    meta_out.to_csv(meta_path, index=False)
    ab_out.to_csv(rel_path, index=True)

    n_total = len(complete)
    print("\n--- exclusion (subjects missing >= 1 common timepoint) ---")
    print(f"subjects: {n_total} total -> kept {len(kept)}, excluded {len(excluded)}")
    print(f"samples : {n_total * 7} total -> kept {len(kept) * 7}, excluded {len(excluded) * 7}")
    print(f"excluded subjects ({len(excluded)}): {excluded}")
    print(f"\nwrote metadata  -> {meta_path} ({len(meta_out)} rows)")
    print(f"wrote abundance -> {rel_path} ({ab_out.shape[1]} samples)")


if __name__ == "__main__":
    main()
