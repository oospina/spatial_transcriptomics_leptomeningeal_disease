#!/usr/bin/env bash
#
#!!! This is not officially supported by 10x
#
# =============================================================================
# Job Script
# =============================================================================
#
#PBS -N st_sample2_man_al
#PBS -V
#PBS -l nodes=1:ppn=16
#PBS -l mem=128gb
#PBS -o out_visium_count.txt
#PBS -e err_visium_count.txt

spaceranger count --id=sample2_LMM_B2 \
	--transcriptome=refdata-gex-GRCh38-2020-A \
	--fastqs=first_run_visium_mkfastq/outs/fastq_path,second_run_visium_mkfastq/outs/fastq_path \
	--sample=LMM_B2 \
	--image=HandE_images_spranger_input/sample2_LMM_B2.tif \
	--slide=[INSERT_VISIUM_SLIDE_ID_HERE] \
	--area=[INSERT_VISIUM_SLIDE_AREA_HERE] \
	--loupe-alignment=../data/loupe_manual_alignments_input/sample2_LMM_B2.json \
