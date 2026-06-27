# BONUS — CF-infant gut microbiome cohort

**Source:** https://doi.org/10.1038/s41591-019-0714-x

16S species-level relative-abundance profiles and clinical metadata for infants
with cystic fibrosis (CF). Here the cohort is turned into a benchmark for
**imputing missing microbiome timepoints**: each subject's longitudinal samples
are masked under different missingness mechanisms, and imputation methods are
scored on how well they recover the masked values.

## What was done

This is a **single cohort** (no country/group split), processed in three stages.

**`01.input_data/`** — raw inputs: the species profile (`profile_Species.csv`) and
the supplementary Excel with the metadata (two sheets: patient-level and
sample/per-month).

**`02.preparation/`** — all preparation and analysis scripts:
- **combine metadata** — `combine_metadata.py` merges the two Excel sheets into one
  per-sample table (patient-level fields broadcast across timepoints), aligns it to
  the profile, and maps each sample to the fixed **7-month grid** (months
  3, 4, 5, 6, 8, 10, 12 → timepoints 0–6). Subjects that don't have all 7 common
  timepoints are **excluded** from both metadata and abundance.
- **feature analysis** — `analyze_metadata.py` audits each metadata feature
  (missingness, static vs time-varying, correlation with time) to choose valid
  **MAR** drivers.
- **mask generation** — `MCAR/MAR/MNAR_data_prep.py` produce the missingness masks.

**`03.processed/`** — generated outputs only:
`01.metadata/`, `02.relative_abundance/`, `03.clr/` (CLR-transformed), `04.sorted/`
(DeepMicroGen input layout), and one folder per mechanism (`MCAR/`, `MAR/`, `MNAR/`)
holding the masks at 9 missingness ratios (0.1–0.9) × 10 replicates.

## Missingness mechanisms

| Mechanism | Meaning | How it's generated |
|---|---|---|
| **MCAR** | missing completely at random | uniform random masking of timepoints |
| **MAR** | depends on observed covariates | `pyampute`, driven by meaningful observed covariates (e.g. `weight %`) |
| **MNAR** | depends on the timepoint itself | direct longitudinal-dropout weighting (later visits more likely missing) |

A mask file is a `subject × timepoint` 0/1 matrix (`0` = masked). The masked
abundance is what each competitor method (DeepMicroGen, BAMITA, …) tries to impute;
run scripts live under `03.run_competitors/`.
