#!/bin/bash
#SBATCH --job-name radEmu_t1e_sims
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:59
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1-900
#SBATCH --output t1e_radEmu_%j.out
#SBATCH --error t1e_radEmu_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R GSL MPFR CMake
Rscript radEmu/scripts/t1e_simulations_updated.R $SLURM_ARRAY_TASK_ID