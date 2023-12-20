#!/usr/bin/env bash
#
#!!! This is not officially supported by 10x
#
# =============================================================================
# Job Script
# =============================================================================
#
#PBS -N st_mkfastq
#PBS -V
#PBS -l nodes=1:ppn=16
#PBS -l mem=128gb
#PBS -o out_visium_mkfastq.txt
#PBS -e err_visium_mkfastq.txt

module load bcl2fastq/2.20 

spaceranger mkfastq --id=first_run_visium_mkfastq --run=data/210524_VH00250_68_AAAKFHCM5 --csv=data/spaceranger_visium_mkfastq.csv --qc
