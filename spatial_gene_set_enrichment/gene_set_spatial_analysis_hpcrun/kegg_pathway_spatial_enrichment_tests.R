# Test for spatial KEGG pathway enrichment (using log data)
library('tidyverse')
library('msigdbr')

cmdargs = commandArgs(trailingOnly=TRUE)
cmdargs = unlist(strsplit(cmdargs, ','))

setwd('.')

# Get gene sets for biological processes
gene_sets = msigdbr(species='Homo sapiens')
gene_sets = filter(gene_sets, gs_subcat == "CP:KEGG") %>% 
  split(x=.$gene_symbol, f=.$gs_name)

reps = 1000
min_units = as.integer(cmdargs[1]) # Minimum number of spots with high expression of a pathway for that pathway to be considered
num_sds = as.double(cmdargs[2])
cores =  as.integer(cmdargs[3])
min_genes = 10
tissue = 'all_spots'
out_dir = '.'
seed = 12345
pval_adj_method = 'BH'

# Read data in Seurat format to get image x,y and SCTransform (same Seurat objects used as input for LDA model selection)
seuratobjs = readRDS('[NAMED_LIST_SEURAT_OBJECTS]')

pval_dfs = list()
set.seed(seed)
for(i in names(seuratobjs)){
  cat(crayon::yellow(paste0('\tSample: ', i, '...\n')))
  # Extract spots to be used in analysis
  # This selection implemented proactively as analysis might later be applied to tissue niches within samples
  tissue_spots = rownames(seuratobjs[[i]]@meta.data)
  
  # Extract gene expression
  exp = as.data.frame(Seurat::GetAssayData(seuratobjs[[i]], assay='Spatial', slot='data'))
  exp = exp[, tissue_spots]
  
  # Extract spot coordinates and match order of spots
  coords_df = seuratobjs[[i]]@images[[i]]@coordinates %>%
    rownames_to_column(var='libname') %>%
    select(libname, xpos=imagecol, ypos=imagerow)
  coords_df = coords_df[match(colnames(exp), coords_df[['libname']]), ]
  
  # Perform tests in parallel
  pval_ls = parallel::mclapply(seq(names(gene_sets)), function(pw){
    pw_genes = gene_sets[[pw]]
    pw_genes = pw_genes[pw_genes %in% rownames(exp)] # Filter genes that aren't in the expression matrix.
    
    # Test if genes in data set are enough
    if(length(pw_genes) >= min_genes){
      # Average expression of genes within pathway for each spot
      pw_avg_exp = apply(exp[pw_genes, ], 2, mean)
      
      # Find spots that highly express this pathway (mean + stdev_t*std in this case)
      avg = mean(pw_avg_exp)
      std = sd(pw_avg_exp)
      exp_thresh = avg + (num_sds*std)
      
      # Extract spots with average expression above sd threshold
      high_spots_bc = names(which(pw_avg_exp >= exp_thresh))
      
      # Are there at least X number of cells/spots?
      if(length(high_spots_bc) >= min_units){
        
        # Compute the distance between these spots' XY coordinates with high expression
        coords_high_spots = coords_df[coords_df[['libname']] %in% high_spots_bc, ]
        distances_high_spots = as.matrix(rdist::rdist(coords_high_spots[, c('xpos', 'ypos')], metric='euclidean')) # Distance computation
        distances_high_spots[upper.tri(distances_high_spots)] = 0 # Make upper half zero to avoid sum of distances twice
        sum_high_distances = sum(distances_high_spots)
        
        # Compute random distance permutations
        sum_rand_distances = c()
        for(rep in 1:reps){
          rand_idx = sample(x=1:nrow(coords_df), size=length(high_spots_bc))
          rand_coord = coords_df[rand_idx, ]
          rand_dist = as.matrix(rdist::rdist(rand_coord[, c('xpos', 'ypos')], metric='euclidean')) # Distance computation for random spots
          rand_dist[upper.tri(rand_dist)] = 0 # Make upper half zero to avoid sum of distances twice
          sum_rand_distances = append(sum_rand_distances, sum(rand_dist))
        }
        
        # Compute p-value
        # Ho: sum of observed distances is larger than sum of random distances
        # Ha: sum of observed distances is smaller than sum of random distances
        #count_test = sum(sum_high_distances > sum_rand_distances) # Times observed dist was higher than random dists
        count_test = sum(sum_rand_distances < sum_high_distances)
        p_val = count_test / reps
        
      } else{
        p_val=NA
      } # Close test minimum spots
    } else{
      p_val=NA
    } # Close test minimum genes
    
    # Get results in data frame form
    pval_tmp = tibble::tibble(gene_set=names(gene_sets[pw]),
                              size_test=length(pw_genes),
                              size_gene_set=length(gene_sets[[pw]]),
                              p_value=p_val)
    
    return(pval_tmp)
    
  }, mc.cores=cores) # Close mclappy
  
  # Compile results in a single data frame
  pval_dfs[[i]] = dplyr::bind_rows(pval_ls) %>%
    tibble::add_column(sample_name=i, .before=1)
  
  # Adjust p-values
  pval_dfs[[i]][['adj_p_value']] = p.adjust(pval_dfs[[i]][['p_value']], method=pval_adj_method)
  pval_dfs[[i]] = pval_dfs[[i]][order(pval_dfs[[i]][['adj_p_value']]), ]
  
}
# Write test results
openxlsx::write.xlsx(pval_dfs, paste0(tissue, '_minspots_', as.character(min_units), '_sd_', as.character(num_sds),'_kegg_pathway_tests.xlsx'))

cat(crayon::green("STenrich completed."))

  
