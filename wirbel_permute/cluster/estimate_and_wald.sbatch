#!/bin/bash
#SBATCH --job-name estimate_and_wald
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:00
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-100
#SBATCH --output estimate_and_wald_%j.out
#SBATCH --error estimate_and_wald_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
Rscript radEmu/wirbel_permute/estimate_and_wald.R $SLURM_ARRAY_TASK_ID