#!/bin/bash
#SBATCH --job-name score_left
#SBATCH --partition medium
#SBATCH --ntasks 1
#SBATCH --time 5-11:00
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-608
#SBATCH --output score_%j.out
#SBATCH --error score_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
Rscript radEmu/wirbel_permute/score_left.R $SLURM_ARRAY_TASK_ID