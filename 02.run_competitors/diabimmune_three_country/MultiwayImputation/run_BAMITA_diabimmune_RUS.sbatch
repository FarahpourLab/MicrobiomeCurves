#!/bin/bash -l
#SBATCH -J bamita_diab_RUS
#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=4-00:00:00
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alireza.dostmohammadi@uni-due.de

cd $SLURM_SUBMIT_DIR

source ~/.bashrc

# Make R + the BAMITA packages available (R module or conda env with the
# packages listed at the top of run_bamita_diabimmune.R).
# module load R/4.3.1
eval "$(conda shell.bash hook)"
conda activate bamita

Rscript run_bamita_diabimmune.R RUS
