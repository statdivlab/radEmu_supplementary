#!/bin/bash
#SBATCH --job-name estimate_and_wald33
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:00
#SBATCH --mem-per-cpu=16000
#SBATCH --output estimate_and_wald33_%j.out
#SBATCH --error estimate_and_wald33_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
Rscript radEmu/wirbel_permute/estimate_and_wald33.R 