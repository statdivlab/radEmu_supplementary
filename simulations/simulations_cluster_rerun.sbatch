#!/bin/bash
#SBATCH --job-name radEmu_sims
#SBATCH --partition medium
#SBATCH --ntasks 1
#SBATCH --time 5-11:59
#SBATCH --mem-per-cpu=12000
#SBATCH --array=1-882
#SBATCH --output rerun_sim_radEmu_%j.out
#SBATCH --error rerun_sim_radEmu_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
module load GSL/2.7-GCC-13.2.0
module load MPFR CMake
Rscript radEmu/scripts/simulations_cluster_rerun.R $SLURM_ARRAY_TASK_ID