
cellranger_folder <- "/data/sc-10x/data-runs/171120-scheele-adipose/agg-180504-unnormalized/outs/filtered_gene_bc_matrices_mex/hg19"
metadata_sheet <- "../files/180406-Cell ID and 10x sample Index_final-extracted columns.txt"
vars.to.regress <- c("nUMI", "percent.mito")
n.pcs <- 15
resolutions <- c(0.5, 0.7, 1, 1.5)
output_beforeqc <- '../output/10x-180504-beforeqc'
output <- "../output/10x-180504"

get_metadata_df <- function(seurobj){
  samples_info <- read.table(metadata_sheet, sep='\t', header=T, stringsAsFactors=F)
  samples_info.ordered <- rbind(samples_info[13:14,], samples_info[1:12,])

  #get the metadata
  sample_names <- samples_info.ordered$Sample.name..Single.Cell.ID.
  diff <- samples_info.ordered$Diff
  ucp1.ctrl <- samples_info.ordered$UCP1.ctrl
  ucp1.ne <- samples_info.ordered$UCP1.NE
  bmi <- samples_info.ordered$BMI
  age <- samples_info.ordered$AGE

  #get the sample names consistent
  sample_names[1] <- 'Supra_4'
  sample_names[2] <- 'Subq_4'

  #get depot labeling (supra, subq, peri or visce)
  depots <- unlist(lapply(sample_names, function(x){
    return(substring(x, 0, nchar(x)-2))
  }))

  #extract the indices from the cellnames
  sample_agg_idx <- as.numeric(sapply(strsplit(seurobj@cell.names, split = "-"), '[[', 2))

  #create metadata df
  df.metadata <- data.frame(row.names=seurobj@cell.names,
                            depot=depots[sample_agg_idx],
                            sample_name=sample_names[sample_agg_idx],
                            diff=diff[sample_agg_idx],
                            ucp1.ctrl=ucp1.ctrl[sample_agg_idx],
                            ucp1.ne=ucp1.ne[sample_agg_idx],
                            bmi=bmi[sample_agg_idx],
                            age=age[sample_agg_idx])
  return(df.metadata)
}


get_filtered_dataset <- function(seurobj){
  seurobj <- FilterCells(seurobj, subset.names=c("nGene", "percent.mito", "nUMI"), low.thresholds = c(200, -Inf, -Inf), high.thresholds = c(9000, 0.08, 110000))
  return(seurobj)
}
