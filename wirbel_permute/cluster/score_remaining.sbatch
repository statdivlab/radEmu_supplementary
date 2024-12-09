#!/bin/bash
#SBATCH --job-name score_remain
#SBATCH --partition long
#SBATCH --ntasks 1
#SBATCH --time 14-11:59
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-180
#SBATCH --output remain_score_%j.out
#SBATCH --error remain_score_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=teichs@uw.edu

module load R-bundle-CRAN
Rscript radEmu/wirbel_permute/score_remaining.R $SLURM_ARRAY_TASK_ID