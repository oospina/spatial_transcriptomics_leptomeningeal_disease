#!/bin/bash
#SBATCH --job-name=stgradient
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --time=200:00:00
#SBATCH --output dist2ref_%A_%a.out
#SBATCH --error dist2ref_%A_%a.err
#SBATCH -a 1-96%12

cd .

module load R

args=$(cat ../../data/arguments_runs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Start analysis: "
echo $(date)
echo ${args}

Rscript correlations_dist2ref_stdeconvolve_collapsed.R $args

cd .

echo "End analysis: "
echo $(date)
