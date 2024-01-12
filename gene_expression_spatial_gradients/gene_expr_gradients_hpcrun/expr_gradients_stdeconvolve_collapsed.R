# Correlation of gene expression with distance to a cluster (based on STdeconvolve 
# collapsed domains) using log-transformed counts (HPC script)

# Load libraries
library('tidyverse')
library('Seurat')

# Capture arguments from command line
cmdargs = commandArgs(trailingOnly=TRUE)

# SPECIFY OPTIONS FOR ANALYSIS
# PARSE ARGUMENTS FROM COMMAND LINE
# Samples to include in analysis. All if empty
samplenames = c()
# Top variable genes to test
topvariable = 2000
# Reference cluster. Distances calculated from that cluster
ref = grep('\\-\\-ref', cmdargs, value=T)
if(length(ref) > 0){
  ref = strsplit(ref, split='=')[[1]][[2]]
}
# Cluster to exclude from analysis? Empty if none
exclude = grep('\\-\\-exclude', cmdargs, value=T)
if(length(exclude) > 0){
  exclude = strsplit(exclude, split='=')[[1]][[2]]
} else{
  exclude = c()
}
# Remove outliers? T/F
out_rm = grep('\\-\\-out_rm', cmdargs, value=T)
if(length(out_rm) > 0){
  out_rm = as.logical(strsplit(out_rm, split='=')[[1]][[2]])
}
# Limit of distances to calculate correlations
# To include all spots, leave empty
limit = grep('\\-\\-limit', cmdargs, value=T)
if(length(limit) > 0){
  limit = as.numeric(strsplit(limit, split='=')[[1]][[2]])
} else{
  limit = c()
}
# Define the function to summarize spot distances (minimum or average)
distsumm = grep('\\-\\-distsumm', cmdargs, value=T)
if(grepl('min|avg', distsumm)){
  distsumm = strsplit(distsumm, split='=')[[1]][[2]]
}
# Specify minimum number of neighbors for reference spot to be included in analysis
min_nb = grep('\\-\\-min_nb', cmdargs, value=T)
if(length(min_nb) > 0){
  min_nb = as.integer(strsplit(min_nb, split='=')[[1]][[2]])
}
# Robust regression? T/F
robust = grep('\\-\\-robust', cmdargs, value=T)
if(length(robust) > 0){
  robust = as.logical(strsplit(robust, split='=')[[1]][[2]])
}

# Make sure no outliers are removed if robust regression
if(robust){
  out_rm = F
}

# Read Seurat objects
# Remove excluded samples for downstream processing
seurat_obj = readRDS(file='[NAMED_LIST_SEURAT_OBJECTS]')

# Get sample names in case no specific sample was selected
if(is.null(samplenames)){
  samplenames = names(seurat_obj)
}

# Read spot annotations for tissue regions and add to Seurat metadata
annots_fp = '../../data/spot_domain_majority_classes.xlsx'
for(i in samplenames){
  annots_tmp = readxl::read_excel(annots_fp, sheet=i)
  # Add spot annotations to meta data
  seurat_obj[[i]]@meta.data = seurat_obj[[i]]@meta.data %>%
    rownames_to_column(var='barcode') %>%
    left_join(., annots_tmp %>%
                select(barcode, domain), by='barcode') %>%
    column_to_rownames(var='barcode')
}
rm(annots_fp, annots_tmp) # Clean environment

# Read distance spreadsheets
ds_fp = paste0('../../data/euclidean_distances_all_spots.xlsx')
dists_dfs = lapply(names(seurat_obj), function(i){
  df_tmp = as.data.frame(readxl::read_excel(ds_fp, sheet=i))
  # Column names and rows should have the same order
  rownames(df_tmp) = colnames(df_tmp)
  return(df_tmp)
})
names(dists_dfs) = names(seurat_obj)
rm(ds_fp) # Clean environment

# Save spots in the different categories (ref, nonref, excl)
spot_cats = lapply(names(seurat_obj), function(i){
  ref_tmp = rownames(seurat_obj[[i]]@meta.data)[seurat_obj[[i]]@meta.data[['domain']] == ref]
  excl_tmp = rownames(seurat_obj[[i]]@meta.data)[seurat_obj[[i]]@meta.data[['domain']] == exclude]
  nonref_tmp = rownames(seurat_obj[[i]]@meta.data)[!(seurat_obj[[i]]@meta.data[['domain']] %in% c(ref, exclude))]  
  return(list(
    ref=ref_tmp,
    nonref=nonref_tmp,
    excl=excl_tmp
  ))
})
names(spot_cats) = names(seurat_obj)

# Identify spots to be removed from reference if not enough neighbors
# Remove sample if no clusters of reference spots
sample_rm = c()
nbs_keep = list()
for(i in samplenames){
  # Get minimum distance among all spots within a sample (for Visium would be approximately the same for any sample)
  min_sample = min(as.data.frame(dists_dfs[[i]][lower.tri(dists_dfs[[i]])]))
  # Get distances among reference spots 
  dists_tmp = dists_dfs[[i]][spot_cats[[i]][["ref"]], spot_cats[[i]][['ref']], drop=F]
  
  # Get number of neighbors within minimum distance
  # NOTE: When dealing with other technologies like SMI, will need to be more flexible with 
  # minimum distances as not an array of equally distant spots. In this case, allowed a "buffer"
  # of a quarter of the minimum distance
  nbs = colSums(dists_tmp >= min_sample * 0.75 & dists_tmp <= min_sample * 1.25)
  if(sum(nbs >= min_nb) < 1){ # At least 1 cluster of spots to include samplenin analysis
    sample_rm = append(sample_rm, i)
  } else{
    nbs_keep[[i]] = names(nbs)[nbs >= min_nb] # Save spots to be kept (enough neighbors)
  }
  rm(list=grep("nbs$|dists_tmp|min_sample", ls(), value=T)) # Clean environment
}
# Remove samples that did not have enough neighbor reference spots
if(length(sample_rm) != 0){
  seurat_obj = seurat_obj[!grepl(paste0(sample_rm, collapse='|'), names(seurat_obj))]
  dists_dfs = dists_dfs[!grepl(paste0(sample_rm, collapse='|'), names(dists_dfs))]
  spot_cats = spot_cats[!grepl(paste0(sample_rm, collapse='|'), names(spot_cats))]
  samplenames = samplenames[!grepl(paste0(sample_rm, collapse='|'), samplenames)]
}
rm(sample_rm) # Clean environment

# Get summarized distances from the reference for each spot in the non-reference
sample_rm = c()
for(i in samplenames){
  # Select spots in analysis (non reference in rows, reference in columns)
  df_tmp = as.data.frame(dists_dfs[[i]][spot_cats[[i]][["nonref"]], spot_cats[[i]][['ref']], drop=F])
  # Remove columns corresponding to spots without enough neighbors
  df_tmp = df_tmp[, colnames(df_tmp) %in% nbs_keep[[i]], drop=F]
  
  # Check that distances are available for the comparison
  # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
  if(nrow(df_tmp) > 1 & ncol(df_tmp) > 0){
    if(distsumm == 'avg'){ # Summarize non-reference spots using minimun or mean distance
      dists_tmp = tibble(barcode=rownames(df_tmp),
                         dist2ref=apply(df_tmp, 1, mean))
    } else{
      dists_tmp = tibble(barcode=rownames(df_tmp),
                         dist2ref=apply(df_tmp, 1, min))
    }
    
    # Add distances to Seurat meta data
    seurat_obj[[i]]@meta.data = seurat_obj[[i]]@meta.data %>%
      rownames_to_column(var='barcode') %>%
      left_join(., dists_tmp, by='barcode') %>%
      column_to_rownames(var='barcode')
    
    rm(dists_tmp) # Clean environment)
  } else{
    # Save samples with all NA distances for removal from analysis
    sample_rm = append(sample_rm, i)
  }
  rm(df_tmp) # Clean environment
}
# Remove samples if all distances were NA (probably sample is all reference)
if(length(sample_rm) != 0){
  seurat_obj = seurat_obj[!grepl(paste0(sample_rm, collapse='|'), names(seurat_obj))]
  samplenames = samplenames[!grepl(paste0(sample_rm, collapse='|'), samplenames)]
}
rm(dists_dfs, sample_rm, nbs_keep, spot_cats) # Clean environment

# Remove distances if outside user-specified limit
if(!is.null(limit)){
  for(i in samplenames){
    # Get lower and upper distance limits 
    # If lower limit is higher than user limit, then set lower limit as upper limit
    if(!all(is.na(seurat_obj[[i]]@meta.data[['dist2ref']]))){
      dist2reflower = min(seurat_obj[[i]]@meta.data[['dist2ref']], na.rm=T)
      if(dist2reflower > limit){
        dist2refupper = dist2reflower
      } else{
        dist2refupper = limit
      }
      # Make NA the distances outside range
      seurat_obj[[i]]@meta.data = seurat_obj[[i]]@meta.data %>%
        mutate(dist2ref=case_when(dist2ref <= dist2refupper ~ as.numeric(dist2ref)))
      
      rm(dist2reflower, dist2refupper) # Clean environment
    }
  }
} 

# Get expression from variable genes
# Genes are identified within the range limit
vargenes_expr = list()
sample_rm = c()
for(i in samplenames){
  # Extract expression data (non transformed counts to be passed to FindVariableFeatures)
  DefaultAssay(seurat_obj[[i]]) = 'Spatial'
  df_tmp = as.data.frame(GetAssayData(seurat_obj[[i]], assay='Spatial', slot='counts')) %>%
    rownames_to_column(var='gene_name') %>%
    filter(!stringr::str_detect(gene_name, '^MT\\-|^RP[S|L]')) %>% # Remove mitochondrial and ribosomal genes
    column_to_rownames(var='gene_name') %>%
    as.matrix()
  
  # Get spots that have at least 1 distance value
  # However, if only one spot, then sample will be removed from analysis as cannot detect variable genes from single spot
  df_tmp = df_tmp[, rownames(seurat_obj[[i]]@meta.data)[!is.na(seurat_obj[[i]]@meta.data$dist2ref)], drop=F]

  # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
  # Variable genes in minimum distance range
  vargenes = c() # In case variable genes cant be detected (less than one spot?)
  if(ncol(df_tmp) > 1){
    vargenes = FindVariableFeatures(df_tmp, verbose=F) %>%
      rownames_to_column(var='gene') %>%
      arrange(desc(vst.variance.standardized)) %>%
      select(gene) %>%
      unlist() %>%
      as.vector()
    vargenes = vargenes[1:topvariable] # Get number of genes defined by user
  }
  rm(df_tmp) # Clean environment
  
  # Get transformed gene expression data (will be used for the correlations with distance)
  # Matrices will contain only the non-reference spots (as defined by non-NA value in distance)
  if(length(vargenes) > 0){
    tr_df = as.matrix(GetAssayData(seurat_obj[[i]], assay='Spatial', slot='data'))
    tr_df = tr_df[, rownames(seurat_obj[[i]]@meta.data)[!is.na(seurat_obj[[i]]@meta.data$dist2ref)], drop=F]
    
    vargenes_expr[[i]] = tr_df[rownames(tr_df) %in% c(vargenes), , drop=F] %>% 
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='barcode') %>%
      left_join(seurat_obj[[i]]@meta.data %>%
                  rownames_to_column(var='barcode') %>%
                  select(barcode, dist2ref), ., by='barcode') %>%
      filter(barcode %in% colnames(tr_df)) %>%
      right_join(seurat_obj[[i]]@images[[i]]@coordinates %>%
                  rownames_to_column(var='barcode') %>%
                  select(barcode, imagerow, imagecol), ., by='barcode')
    
    rm(tr_df, vargenes) # Clean environment
  } else{
    sample_rm = append(sample_rm, i)
  }
}
# Remove samples if all distances were NA (probably sample is all reference)
if(length(sample_rm) != 0){
  seurat_obj = seurat_obj[!grepl(paste0(sample_rm, collapse='|'), names(seurat_obj))]
  samplenames = samplenames[!grepl(paste0(sample_rm, collapse='|'), samplenames)]
}
rm(sample_rm) # Clean environment

# Detect gene expression outlier spots for each sample and gene
if(out_rm & !robust){
  outs_dist2ref = list()
  
  for(i in samplenames){
    outs_dist2ref[[i]] = list()
    
    dfdist2ref = vargenes_expr[[i]][!is.na(vargenes_expr[[i]]$dist2ref), ] %>%
      dplyr::select(-barcode, -imagerow, -imagecol, -dist2ref)

    for(gene in colnames(dfdist2ref)){
      # Calculate gene expression quartiles
      quarts = quantile(dfdist2ref[[gene]], probs=c(0.25, 0.75))
      # Calculate inter-quartile range
      iqr_dist2ref = IQR(dfdist2ref[[gene]])
      # Calculate distribution lower and upper limits
      low_up_limits = c((quarts[1]-1.5*iqr_dist2ref), 
                        (quarts[2]+1.5*iqr_dist2ref))
      
      # Save outliers (barcodes)
      outs_dist2ref[[i]][[gene]] = dfdist2ref[['barcode']][ dfdist2ref[[gene]] < low_up_limits[1] | 
                                                              dfdist2ref[[gene]] > low_up_limits[2] ]
    }
  }
  rm(list=grep("iqr|quarts|low_up|dfdist2ref", ls(), value=T)) # Clean environment
}

# Calculate Spearman's correlations
dist_cor = list()
for(i in samplenames){
  # Initialize data frame to store results
  dist_cor[[i]] = tibble(gene=character(),
                         dist2ref_est=numeric(),
                         dist2ref_p=numeric(),
                         dist2ref_spearman=numeric(),
                         dist2ref_spearman_p=numeric())
  
  # CORRELATIONS DISTANCE TO REFERENCE CLUSTER
  genes_sample = colnames(vargenes_expr[[i]] %>% dplyr::select(-barcode, -imagerow, -imagecol, -dist2ref))
  for(gene in genes_sample){
    df_gene = vargenes_expr[[i]] %>%
      dplyr::select(barcode, dist2ref, .data[[gene]])
    
    lm_res = list(estimate=NA, estimate_p=NA) 
    cor_res = list(estimate=NA, p.value=NA) 
    if(out_rm & !robust){ # Regular linear models after removal of outliers
      # Remove outliers
      df_gene_outrm = df_gene %>% 
        filter(!(barcode %in% outs_dist2ref[[i]][[gene]]))
      if(nrow(df_gene_outrm) > 1){
        # Run linear model and get summary
        lm_tmp = lm(df_gene_outrm[[gene]] ~ df_gene_outrm[['dist2ref']])
        lm_summ_tmp = summary(lm_tmp)[['coefficients']]
        if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
          lm_res = list(estimate=lm_summ_tmp[2,1],
                        estimate_p=lm_summ_tmp[2,4])
        }
        # Calculate Spearman correlation
        cor_res = cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman')
      }
    } else {
      if(robust){ # Robust linear models?
        df_gene_range = df_gene
        if(nrow(df_gene_range) > 1){
          # Run robust linear model and get summary
          lm_tmp = MASS::rlm(df_gene_range[[gene]] ~ df_gene_range[['dist2ref']], maxit=100)
          if(lm_tmp[['converged']] & lm_tmp[['coefficients']][2] != 0){ # Check the model converged and an effect was estimated
            # Run Wald test (MASS::rlm does not provide a p-value)
            lm_test_tmp = sfsmisc::f.robftest(lm_tmp)
            lm_res = list(estimate=summary(lm_tmp)[['coefficients']][2,1],
                          estimate_p=lm_test_tmp[['p.value']])
            # Calculate Spearman correlation
            cor_res= cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')
          }
        }
      } else{ # Regular linear models without outlier removal
        df_gene_range = df_gene
        if(nrow(df_gene_range) > 1){
          lm_tmp = lm(df_gene_range[[gene]] ~ df_gene_range[['dist2ref']])
          lm_summ_tmp = summary(lm_tmp)[['coefficients']]
          if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
            lm_res = list(estimate=lm_summ_tmp[2,1],
                          estimate_p=lm_summ_tmp[2,4])
          }
          cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')
        }
      }
    }
    
    # Create row with results
    tibble_tmp = tibble(gene=gene, 
                        dist2ref_est=lm_res$estimate,
                        dist2ref_p=lm_res$estimate_p,
                        dist2ref_spearman=cor_res$estimate,
                        dist2ref_spearman_p=cor_res$p.value)
    
    rm(list=grep("lm_|_res|_test|cor_|df_gene", ls(), value=T)) # Clean environment
    
    # Add row to result table if there is at least one correlation value
    if(sum(!is.na(tibble_tmp)) > 1){
      dist_cor[[i]] =  bind_rows(dist_cor[[i]], tibble_tmp)
      rm(tibble_tmp) # Clean environment
    }
  }
  
  rm(genes_sample) # Clean environment
  
  # Adjust p-values for multiple comparison
  dist_cor[[i]] = dist_cor[[i]] %>%
    add_column(adjpval_dist2ref=p.adjust(.$dist2ref_spearman_p)) %>%
    relocate(adjpval_dist2ref, .after=dist2ref_spearman_p) %>%
    arrange(adjpval_dist2ref)
}

# Save correlations to file
fname = paste0('correlations_geneexpr_with_dist_ref_', 
               paste0(ref, collapse='_'), '_', distsumm, '_outrm_', out_rm)
if(!is.null(limit)){
  fname = paste0(fname, '_limitrange_', limit)
}
if(!is.null(exclude)){
  fname = paste0(fname, '_excl_', paste0(exclude))
}
if(!is.null(min_nb)){
  fname = paste0(fname, '_minnb_', paste0(min_nb))
}
if(!is.null(robust)){
  fname = paste0(fname, '_robreg_', paste0(robust))
}
fname = paste0(fname, '.xlsx')
openxlsx::write.xlsx(dist_cor, file=fname)

