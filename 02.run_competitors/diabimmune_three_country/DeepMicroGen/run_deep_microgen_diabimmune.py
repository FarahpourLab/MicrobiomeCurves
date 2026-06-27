"""DeepMicroGen runner for the DIABIMMUNE "three country" dataset.

One COUNTRY per job (passed as ``argv[1]``: EST / FIN / RUS). NO hyperparameter
tuning -- a single fixed learning rate / dropout / epochs. For that country, and
for every missingness mechanism (MCAR, MAR, MNAR), ratio and replicate, it calls
``deepMicroGen.py`` on the sorted CLR profile + the mask so the model imputes the
masked timepoints.

Input  : 01.Benchmarking Datasets/DIABIMMUNE three country/03.processed/
           04.sorted/sorted_clr_transformation_relative_abundance_l7_species_<C>_<pseudo>.csv
           <MECH>/mask/l7_species/<ratio>/masked_metadata_l7_species_<C>.rep_*.csv
Output : 04.Benchmarking/DIABIMMUNE three country/DeepMicroGen/<MECH>/<ratio>/
"""
import os
import sys
import glob
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATASET = "DIABIMMUNE three country"
TAX_LEVEL = "l7_species"
MECHANISMS = ["MCAR", "MAR", "MNAR"]

# fixed parameters (no tuning)
LEARNING_RATE = "0.001"
DROPOUT_RATE = "0.7"
EPOCHS = "3000"


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python run_deep_microgen_diabimmune.py <COUNTRY>")
    country = sys.argv[1]

    processed = os.path.join(ROOT, "01.Benchmarking Datasets", DATASET, "03.processed")
    sorted_dir = os.path.join(processed, "04.sorted")
    out_base = os.path.join(ROOT, "04.Benchmarking", DATASET, "DeepMicroGen")

    input_files = glob.glob(os.path.join(
        sorted_dir,
        f"sorted_clr_transformation_relative_abundance_{TAX_LEVEL}_{country}_*.csv"))
    if not input_files:
        raise SystemExit(f"no sorted input found for country '{country}' in {sorted_dir}")

    for input_file_path in input_files:
        profiles_name = os.path.basename(input_file_path).split(".csv")[0]

        for mechanism in MECHANISMS:
            mask_glob = os.path.join(
                processed, mechanism, "mask", TAX_LEVEL, "*",
                f"masked_metadata_{TAX_LEVEL}_{country}.rep_*.csv")

            for mask_path in sorted(glob.glob(mask_glob)):
                ratio = os.path.basename(os.path.dirname(mask_path))
                out_dir = os.path.join(out_base, mechanism, ratio)
                os.makedirs(out_dir, exist_ok=True)

                mask_name = ".".join(os.path.basename(mask_path).split(".csv")[0].split(".")[1:])
                out_file = os.path.join(out_dir, f"imputed_{profiles_name}_{mask_name}.csv")
                if os.path.isfile(out_file):
                    continue

                cmd = ["python", "deepMicroGen.py",
                       input_file_path, mask_path, out_dir,
                       LEARNING_RATE, DROPOUT_RATE, EPOCHS]
                subprocess.run(cmd, check=True)
