#!/bin/bash -l
#SBATCH -J dmg_diab_FIN
#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=4-00:00:00
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alireza.dostmohammadi@uni-due.de

cd $SLURM_SUBMIT_DIR

source ~/.bashrc

# Proper conda activation
eval "$(conda shell.bash hook)"
conda activate deep_microgen

python run_deep_microgen_diabimmune.py FIN
