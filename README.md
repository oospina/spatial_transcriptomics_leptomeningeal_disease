# Analysis of spatial transcriptomics data from melanoma leptomeningeal (LMM) and extra-cranial metastases

This repository contains scripts associated with the pre-print manuscript entitled: *Spatial 
transcriptomics analysis identifies a unique tumor-promoting function of the meningeal stroma in 
melanoma leptomeningeal disease*. The manuscript pre-print can be accessed at [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.18.572266v1)
or [CellPress Sneak Peek](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4685391).
Data used in this code is stored in the Gene Expression Omnibus (GEO) (accession number available
upon publication). Whenever possible and/or relevant, intermediate or result files are 
included in this repository. Unix and R code was contributed by Oscar Ospina.

The .fastq files from Illumina outputs were processed using [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/tutorials/count-ff-tutorial). 
The [Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html) package was used to 
import and further process gene expression counts. Biological identification of the spots was 
achieved using [STdeconvolve](https://jef.works/STdeconvolve/). Gene expression gradients (STgradient) and 
spatial gene set enrichment analysis (STenrich) were conducted in [spatialGE](https://fridleylab.github.io/spatialGE/articles/spatial_enrichment_gradients_smi.html).

## `pre_processing`
The `pre_processing` folder contains scripts to use `spaceranger mkfastq` in order to 
generate .fastq files from the Illumina outputs and `spaceranger count` to generate gene
expression count matrices.
* `pre_processing/visium_mkfastq.sh`: A PBS script to run `spaceranger mkfastq` on an
HPC environment and obtain .fastq files from the Illumina sequencing runs. The script was 
executed separately for each Illumina run. The resulting folder contains demultiplexed .fastq 
files for each sample.
* `pre_processing/visium_count.sh`: An example PBS script to run `spaceranger count` on an
HPC environment and obtain gene expression counts. Counts are assigned to genes in the GRCh38-2020-A
reference genome (downloadable at https://www.10xgenomics.com/support/software/cell-ranger/downloads/).
In order to run the script, a folder containing the tissue images ("HandE\_images\_spranger\_input") 
taken on each Visium slide area.

## `unsupervised_clustering`
The `unsupervised_clustering` folder contains code to perform unsupervised cluster detection and
differential gene expression of the samples (separately and integrated).
* `unsupervised_clustering/unsuperv_cluster_merged.Rmd`: A pipeline to read Space Ranger outputs and 
create a Seurat object with gene expression data from all the samples merged in a single 
matrix. Both log-transformed and SCT-transformed data are generated. Louvain clusters are
identified.
* `unsupervised_clustering/unsuperv_cluster_separate.Rmd`: A pipeline to read Space Ranger outputs and 
create a Seurat object for each sample, resulting in a named list of seurat objects. Both 
log-transformed and SCT-transformed data are generated. Louvain clusters are identified for each sample.
The named list containing Seurat objects is used in downstream analyses (e.g., spot phenoytyping), and
is referenced as '[NAMED_LIST_SEURAT_OBJECTS]' in the scripts.

## `spot_phenotyping_stdeconvolve`
The `spot_phenotyping_stdeconvolve` folder contains code to fit Latent Dirichlet Allocation (LDA)
models using STdeconvolve. The pipeline involves fitting a series models with varying number of topics,
selecting the bets-fitting model, and calculating GSEA scores to assign biological identities to 
the topics (based on the BluePrint data base). The scripts are run in the following order:
* `spot_phenotyping_stdeconvolve/stdeconvolve_model_fitting.Rmd`: Code to fit a series of LDA
models, one sample at a time.
* `spot_phenotyping_stdeconvolve/stdeconvolve_selection_k.Rmd`: Code to select the best-fitting
model for each sample (i.e., find the most likely number of topics in the sample). The procedure is
done one sample at a time.
* `spot_phenotyping_stdeconvolve/gsea_topic_deconvolution.Rmd`: Code to automatically assign 
biological identities to the topics using GSEA. The results are manually curated after inspecting
log-fold changes of key genes in each topic (see `manuscript_figures`).

## `spatial_gene_set_enrichment`
The `spatial_gene_set_enrichment` folder contains code to conduct spatial gene set enrichment
using the STenrich algorithm. Details about the STenrich algorithm can be found [here](https://fridleylab.github.io/spatialGE/articles/spatial_enrichment_gradients_smi.html).
The code was written to be used in a high performance computing (HPC) environment, however, 
interested readers can also run the STenrich algorithm using the function implemented in spatialGE.
* `spatial_gene_set_enrichment/gene_set_spatial_analysis_hpcrun/run_kegg_pathway_spatial_enrichment_tests_minunits_20_numsds_1.sh`:
A Slurm script to submit a run of the `kegg_pathway_spatial_enrichment_tests.R` script to an
HPC environment. A similar script is provided to run STenrich using GO gene sets. Note: The 
STenrich algorithm can also be executed on a regular computer, using [spatialGE](https://fridleylab.github.io/spatialGE/articles/spatial_enrichment_gradients_smi.html).
* `spatial_gene_set_enrichment/gene_set_spatial_analysis_hpcrun/kegg_pathway_spatial_enrichment_tests.R`
An R script that executes the STenrich algorithm. A similar script is provided to run STenrich 
using GO gene sets. Note: The STenrich algorithm can also be executed on a regular computer, 
using [spatialGE](https://fridleylab.github.io/spatialGE/articles/spatial_enrichment_gradients_smi.html).
* `spatial_gene_set_enrichment/spatial_plot_geneset_avgexpr.Rmd`: An R script to generate plots
depicting the average expression of a gene set at each Visium spot. 
* `spatial_gene_set_enrichment/summarize_stenrich_results.Rmd`: Code to generate a summary 
(heatmap) of the gene sets with evidence of significant spatial enrichment.

## `manuscript_figures`
The `manuscript_figures` folder contains a script to generate figures included in the manuscript
* `manuscript_figures/figures_lmm_manuscript.Rmd`: An R Markdown script to generate the individual
plots in the panel figures of the manuscript.

## `data`
The `data` folder contains some of the files necessary to run the spatial analysis pipeline.
Other files need to be produced by the user given constraints in size.
* `data/loupe_manual_alignments_input`: A series of .json files containing coordinates for
manual alignment of spots to tissue images. The manual alignment was done only for samples 
when it was necessary. The files are used during execution of `spaceranger count`.
* `data/deconv_matrices`: Matrices containing relative abundance of each LDA topic at each Visium
spot
* `data/spaceranger_visium_mkfastq.csv`: The .csv file used in `spaceranger mkfastq`,
which that relates sample names and Illumina dual-indexes used in library construction.
* `data/topic_manual_deconvolution.xlsx`: The final biological annotations given to each LDA 
topic and the original annotation given by GSEA. The rationale for changes in final annotation
is also provided in this file.
