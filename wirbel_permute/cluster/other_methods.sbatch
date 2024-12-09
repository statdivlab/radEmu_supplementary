#!/bin/bash
#SBATCH --job-name other_methods
#SBATCH --partition medium
#SBATCH --ntasks 1
#SBATCH --time 3-11:00
#SBATCH --mem-per-cpu=26000
#SBATCH --array=1-100
#SBATCH --output other_methods_%j.out
#SBATCH --error other_methods_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN GSL MPFR CMake
Rscript radEmu/wirbel_permute/other_methods.R $SLURM_ARRAY_TASK_ID