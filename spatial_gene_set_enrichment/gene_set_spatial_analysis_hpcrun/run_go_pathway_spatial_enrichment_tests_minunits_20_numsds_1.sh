#!/bin/bash
#SBATCH --job-name=spatial_go
#SBATCH --cpus-per-task=31
#SBATCH --mem-per-cpu=12G
#SBATCH --time=200:00:00
#SBATCH --output output_go_minunits_20_numsds_1.out
#SBATCH --error errors_go_minunits_20_numsds_1.err

cd .

module load R

echo "Start analysis: "
echo $(date)

# Arguments after script name:
# Minimum high expr spots for pathway to be tested (min_units)
# Number of stdevs to set high expression threshold (num_sds)
# Number of cores to use (cores)
Rscript go_pathway_spatial_enrichment_tests.R 20 1 30

cd .

echo "End analysis: "
echo $(date)
