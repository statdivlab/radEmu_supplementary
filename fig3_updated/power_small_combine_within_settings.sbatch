#!/bin/bash
#SBATCH --job-name combine_power_res
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:59
#SBATCH --mem-per-cpu=15000
#SBATCH --output combine_power_%j.out
#SBATCH --error combine_power_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
Rscript radEmu/scripts/power_small_combine_within_settings.R 