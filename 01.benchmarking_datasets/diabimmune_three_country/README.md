# DIABIMMUNE — three-country cohort

**Source:** https://diabimmune.broadinstitute.org/diabimmune/three-country-cohort

16S species-level relative-abundance profiles and per-sample metadata for infants
from Finland (FIN), Estonia (EST) and Russian Karelia (RUS). Here the cohort is
turned into a benchmark for **imputing missing microbiome timepoints**: each
subject's longitudinal samples are masked under different missingness mechanisms,
and imputation methods are scored on how well they recover the masked values.

## What was done

The data is split **by country** and processed in three stages.

**`01.input_data/`** — raw inputs: the species table (`rel-species-table.csv`) and
the sample metadata (`metadata.csv`).

**`02.preparation/`** — all preparation and analysis scripts:
- **common timepoints** — each country is sampled at slightly different months, so
  for each country we keep the months shared by most subjects (5 common timepoints)
  and drop subjects that don't have all of them, or that are missing a MAR-driver
  covariate. This yields one shared, complete metadata + abundance per country.
- **CLR + sort** — the relative abundance is CLR-transformed and reordered into the
  layout DeepMicroGen expects (OTUs grouped into phylum "clusters").
- **feature analysis** — `analyze_metadata.py` audits each metadata feature
  (missingness, static vs time-varying, correlation with time) to choose valid
  **MAR** drivers.
- **mask generation** — `MCAR/MAR/MNAR_data_prep.py` produce the missingness masks.

**`03.processed/`** — generated outputs only:
`01.metadata/`, `02.relative_abundance/`, `03.clr/`, `04.sorted/`, and one folder per
mechanism (`MCAR/`, `MAR/`, `MNAR/`) holding the masks at 9 missingness ratios
(0.1–0.9) × 10 replicates.

## Missingness mechanisms

| Mechanism | Meaning | How it's generated |
|---|---|---|
| **MCAR** | missing completely at random | uniform random masking of timepoints |
| **MAR** | depends on observed covariates | `pyampute`, driven by a per-country combination of meaningful, non-time-correlated covariates |
| **MNAR** | depends on the timepoint itself | direct longitudinal-dropout weighting (later visits more likely missing) |

A mask file is a `subject × timepoint` 0/1 matrix (`0` = masked). The masked
abundance is what each competitor method (DeepMicroGen, BAMITA, …) tries to impute;
run scripts live under `03.run_competitors/`.
