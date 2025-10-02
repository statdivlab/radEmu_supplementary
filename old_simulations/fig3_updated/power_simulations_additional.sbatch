#!/bin/bash
#SBATCH --job-name radEmu_t1e_sims
#SBATCH --partition medium
#SBATCH --ntasks 1
#SBATCH --time 5-11:59
#SBATCH --mem-per-cpu=8000
#SBATCH --array=1-540
#SBATCH --output power_radEmu_%j.out
#SBATCH --error power_radEmu_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN GSL MPFR CMake
Rscript radEmu/scripts/power_simulations_additional.R $SLURM_ARRAY_TASK_ID