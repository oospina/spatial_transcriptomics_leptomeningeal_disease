# Analysis of spatial transcriptomics data from melanoma leptomeningeal (LMM) and extra-cranial metastases

This repository contains scripts associated with the pre-print manuscript entitled: *Spatial 
transcriptomics analysis identifies a unique tumor-promoting function of the meningeal stroma in 
melanoma leptomeningeal disease*. The manuscript can be accessed [here](https://www.biorxiv.org/content/10.1101/2023.12.18.572266v1). 
Data used in this code is stored in the Gene Expression Omnibus (GEO) under accession 
number (XXXXXXXXX). Whenever possible and/or relevant, intermediate or result files are 
included in this repository. Unix and R code was contributed by Oscar Ospina.

The .fastq files from Illumina outputs were processed using [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/tutorials/count-ff-tutorial). 
The [`Seurat`](https://satijalab.org/seurat/articles/spatial_vignette.html) package was used to 
import and further process gene expression counts. Biological identification of the spots was 
achieved using [`STdeconvolve`](https://jef.works/STdeconvolve/). Gene expression gradients (STgradient) and 
spatial gene set enrichment analysis (STenrich) were conducted in [`spatialGE`](https://fridleylab.github.io/spatialGE/articles/spatial_enrichment_gradients_smi.html)

## `pre_processing`:
The `pre_processing` folder contains scripts to use `spaceranger mkfastq` in order to 
generate .fastq files from the Illumina outputs and `spaceranger count` to generate gene
expression count matrices.
* `pre_processing/code/visium_mkfastq.sh`: A PBS script to run `spaceranger mkfastq` on an
HPC environment and obtain .fastq files from the Illumina sequencing runs. The script was 
executed separately for each Illumina run. The resulting folder contains demultiplexed .fastq 
files for each sample.
* `pre_processing/code/visium_count.sh`: An example PBS script to run `spaceranger count` on an
HPC environment and obtain gene expression counts. Counts are assigned to genes in the GRCh38-2020-A
reference genome (downloadable at https://www.10xgenomics.com/support/software/cell-ranger/downloads/).
In order to run the script, a folder containing the tissue images ("HandE\_images\_spranger\_input") 
taken on each Visium slide area.
* `pre_processing/data/spaceranger_visium_mkfastq.csv`: The .csv file used in `spaceranger mkfastq`,
which that relates sample names and Illumina dual-indexes used in library construction.
* `pre_processing/data/loupe_manual_alignments_input`: A series of .json files containing coordinates for
manual alignment of spots to tissue images. The manual alignment was done only for samples where it
was necessary. The files are used during execution of `spaceranger count`.

