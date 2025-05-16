#!/bin/bash
#SBATCH --job-name combine_t1e_res
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:59
#SBATCH --mem-per-cpu=15000
#SBATCH --output combine_t1e_%j.out
#SBATCH --error combine_t1e_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
Rscript radEmu/scripts/combine_within_settings_linda.R 