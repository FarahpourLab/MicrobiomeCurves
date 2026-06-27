"""DeepMicroGen runner for the BONUS (CF infants) dataset.

Single cohort (one sorted input, no group token), 7 timepoints. NO hyperparameter
tuning -- a single fixed learning rate / dropout / epochs. For every missingness
mechanism (MCAR, MAR, MNAR), ratio and replicate it calls ``deepMicroGen.py`` on
the sorted CLR profile + the mask, and the model imputes the masked timepoints.

Input  : 01.Benchmarking Datasets/BONUS/03.processed/
           04.sorted/sorted_clr_transformation_relative_abundance_l7_species_<pseudo>.csv
           <MECH>/mask/l7_species/<ratio>/masked_metadata_l7_species.rep_*.csv
Output : 04.Benchmarking/BONUS/DeepMicroGen/<MECH>/<ratio>/
"""
import os
import glob
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATASET = "BONUS"
TAX_LEVEL = "l7_species"
MECHANISMS = ["MCAR", "MAR", "MNAR"]

# fixed parameters (no tuning)
LEARNING_RATE = "0.001"
DROPOUT_RATE = "0.7"
EPOCHS = "3000"


if __name__ == "__main__":
    processed = os.path.join(ROOT, "01.Benchmarking Datasets", DATASET, "03.processed")
    sorted_dir = os.path.join(processed, "04.sorted")
    out_base = os.path.join(ROOT, "04.Benchmarking", DATASET, "DeepMicroGen")

    input_files = glob.glob(os.path.join(
        sorted_dir, f"sorted_clr_transformation_relative_abundance_{TAX_LEVEL}_*.csv"))
    if not input_files:
        raise SystemExit(f"no sorted input found in {sorted_dir}")

    for input_file_path in input_files:
        profiles_name = os.path.basename(input_file_path).split(".csv")[0]

        for mechanism in MECHANISMS:
            mask_glob = os.path.join(
                processed, mechanism, "mask", TAX_LEVEL, "*",
                f"masked_metadata_{TAX_LEVEL}.rep_*.csv")

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
