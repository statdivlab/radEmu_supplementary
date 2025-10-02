#!/bin/bash
#SBATCH --job-name radEmu_power_sims
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1-900
#SBATCH --output sm_power_radEmu_%j.out
#SBATCH --error sm_power_radEmu_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN GSL MPFR CMake
Rscript radEmu/scripts/power_simulations_small_updated.R $SLURM_ARRAY_TASK_ID