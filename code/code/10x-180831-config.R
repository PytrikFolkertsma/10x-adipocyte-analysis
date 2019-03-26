cellranger_folder <- "/data/sc-10x/data-runs/171120-scheele-adipose/agg-180831-unnormalized/outs/filtered_gene_bc_matrices_mex/hg19"
vars.to.regress <- NULL
n.pcs <- 21
resolutions <- c(0.5, 0.7, 1, 1.5)
output_beforeqc <- '../output/10x-180831-beforeqc'
output <- '../output/10x-180831'

get_metadata_df <- function(seurobj){
  timepoints <- c('T1', 'T2', 'T3', 'T4', 'T5', 'T6')
  time_combined <- c(1,1,1,2,2)
  #extract the indices from the cellnames
  sample_agg_idx <- as.numeric(sapply(strsplit(seurobj@cell.names, split = "-"), '[[', 2))
  #assign the correct sample names
  df.metadata <- data.frame(row.names=seurobj@cell.names,
                            timepoint=timepoints[sample_agg_idx],
                            time_combined=time_combined[sample_agg_idx]
  )
  return(df.metadata)
}


get_filtered_dataset <- function(seurobj){
  #No filtering for this dataset since we can filter out doublets with
  #Demuxlet, and cells with a high percent.mito might be brown adipocytes.
  #Remove out bad quality sample 6.
  seurobj <- SubsetData(seurobj, cells.use=rownames(seurobj@meta.data)[seurobj@meta.data$timepoint != 'T6'])
  return(seurobj)
}
