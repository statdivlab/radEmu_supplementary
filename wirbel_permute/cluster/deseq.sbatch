#!/bin/bash
#SBATCH --job-name deseq
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:00
#SBATCH --mem-per-cpu=26000
#SBATCH --array=1-100
#SBATCH --output deseq_%j.out
#SBATCH --error deseq_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN GSL MPFR CMake
Rscript radEmu/wirbel_permute/deseq.R $SLURM_ARRAY_TASK_ID