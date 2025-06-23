#!/bin/bash
#SBATCH --job-name ifaa
#SBATCH --partition medium
#SBATCH --ntasks 1
#SBATCH --time 3-11:00
#SBATCH --mem-per-cpu=26000
#SBATCH --array=1-20
#SBATCH --output ifaa_%j.out
#SBATCH --error ifaa_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN GSL MPFR CMake
Rscript radEmu/wirbel_permute/ifaa.R $SLURM_ARRAY_TASK_ID