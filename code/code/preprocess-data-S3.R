
################################################################################
#
# Creates a preprocessed Seurat object of CellRanger output.
# Executes the following steps:
# - Adds all metadata
# - Filters, normalizes and scales data
# - Computes clusters of several resolutions
# - Calculates cell cycle scores and adds them to metadata
# - Performs tSNE
# If a Seurat object is given as input, only clustering, PCA and TSNE are performed.
#
################################################################################

#Parse command line arguments
library(optparse)

option_list <- list(
  make_option(c('-c', '--config'), type='character', help='Path to R file containing dataset specific variables such as metadata, QC cutoffs, num PCs, etc.'),
  make_option(c('-f', '--file'), type='character', help='OPTIONAL. If a path to a Seurat object is provided, only clustering, PCA and tSNE are applied.')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if (is.null(opt$config)){
  print_help(optionparser)
  quit()
}

################################################################################

library(Seurat)
library(magrittr)
library(dplyr)

################################################################################

#Set variables
source(opt$config) #contains dataset specific variables such as cutoffs, n.pcs, etc.
cellcycle_genes <- "../files/regev_lab_cell_cycle_genes.txt"

###############################################################################

#if a Seurat object is given as input, skip QC and only do PCA, clustering and tSNE.
if(is.null(opt$file)){

  print('>>>LOADING DATA')

  df.10x <- Read10X(cellranger_folder)
  seurobj <- CreateSeuratObject(df.10x, project='10x-adipocyte', min.cells = 3, min.features = 200)
  df.metadata <- get_metadata_df(seurobj) #function from config file
  seurobj <- AddMetaData(seurobj, df.metadata)

  print('>>>QUALITY CONTROL')

  print('Calculating and percent.mito...')
  mito.genes <- grep(pattern = "^MT-", x = rownames(seurobj), value = TRUE, ignore.case=TRUE)
  percent.mito <- Matrix::colSums(GetAssayData(object = seurobj, slot = "counts")[mito.genes, ])/Matrix::colSums(GetAssayData(object = seurobj, slot = "counts"))
  seurobj <- AddMetaData(seurobj, metadata=percent.mito, col.name="percent.mito")

  print(paste('Saving object before QC:', output_beforeqc))
  saveRDS(seurobj, output_beforeqc)

  print('Filtering cells...')
  seurobj <- get_filtered_dataset(seurobj) #function from config file

  print('Normalizing data...')
  seurobj <- NormalizeData(seurobj, normalization.method = "LogNormalize", scale.factor = 1e4)

  print('Scaling data...')
  if (is.null(vars.to.regress)){
    seurobj <- ScaleData(seurobj, do.par=T, num.cores=10)
  } else {
    seurobj <- ScaleData(seurobj, vars.to.regress = vars.to.regress, do.par=T, num.cores=10)
  }
  force_snn_recalc = F
} else {
  seurobj <- readRDS(opt$file)
  force_snn_recalc = T
}

print('>>>CLUSTERING')
print('Finding variable genes...')
seurobj <- FindVariableFeatures(seurobj)
print(paste('Nr. of variable genes:', length(VariableFeatures(object = seurobj))))
print('Running PCA...')
seurobj <- RunPCA(seurobj, npcs=50)
print('Running clustering...')

seurobj <- FindNeighbors(seurobj, dims=1:n.pcs)

for (res in resolutions){
  seurobj <- FindClusters(seurobj, reduction.type = "pca", resolution = res)
}


print('>>>CALCULATING CC SCORES')
cc.genes <- readLines(cellcycle_genes)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
seurobj <- CellCycleScoring(seurobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

print('>>>RUNNING TSNE')
seurobj <- RunTSNE(seurobj, reduction.use='pca', dims=1:n.pcs)

print('>>>SAVING DATASET')
if(is.null(opt$file)){
  print(paste('Saving dataset:', output))
  saveRDS(seurobj, output)
} else {
  print(paste('Saving dataset:', opt$file))
  saveRDS(seurobj, opt$file)
}

