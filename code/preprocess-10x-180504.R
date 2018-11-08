
################################################################################
#
# Creates a Seurat object of the 10x-180504 dataset Executes the following steps:
# - Adds all metadata
# - Filters, normalizes and scales data (filters genes 200-9000,
#   percent.mito < 0.08, nUMI < 110000, regresses out nUMI and percent.mito)
# - Computes clusters (using 15 PC's based on elbow plot) of several resolutions
#   and adds them to metadata
# - Calculates cell cycle scores and adds them to metadata
# - Performs tSNE (15 PC's (determined before from elbow plot), default perplexity)
#
################################################################################

#Load required libraries

library(Seurat)
library(magrittr)
library(dplyr)
library(optparse)

################################################################################

#Set variables

cellranger_folder <- "/data/sc-10x/data-runs/171120-scheele-adipose/agg-180504-unnormalized/outs/filtered_gene_bc_matrices_mex/hg19"
metadata_sheet <- "../files/180406-Cell ID and 10x sample Index_final-extracted columns.txt"
cellcycle_genes <- "../files/regev_lab_cell_cycle_genes.txt"
output <- "../output/10x-180504"

n.pcs <- 15
resolutions <- c(0.5, 0.7, 1, 1.5)

################################################################################

print('>>>LOADING DATA')

df.10x <- Read10X(cellranger_folder)
seurobj <- CreateSeuratObject(df.10x, min.cells = 3, min.genes = 200, is.expr = 0)

################################################################################

print('>>>ADDING METADATA')

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

seurobj <- AddMetaData(seurobj, df.metadata)

print('Nr. of cells per sample before QC: ')
print(seurobj@meta.data %>% count(sample_name))

################################################################################

print('>>>QUALITY CONTROL')

print('Calculating percent.mito...')
mito.genes <- grep(pattern = "^MT-", x = rownames(seurobj@data), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(seurobj@raw.data[mito.genes, ])/Matrix::colSums(seurobj@raw.data)
seurobj <- AddMetaData(seurobj, metadata=percent.mito, col.name="percent.mito")

print('Saving object before QC')
saveRDS(seurobj, '../output/10x-180504-beforeQC')

print('Filtering cells...')
seurobj <- FilterCells(seurobj, subset.names=c("nGene", "percent.mito", "nUMI"), low.thresholds = c(200, -Inf, -Inf), high.thresholds = c(9000, 0.08, 110000))

print('Nr. of cells per sample after QC:')
print(seurobj@meta.data %>% count(sample_name))

print('Normalizing data...')
seurobj <- NormalizeData(seurobj, normalization.method = "LogNormalize", scale.factor = 1e4)

print('Scaling data (regressing out nUMI and percent.mito)...')
seurobj <- ScaleData(seurobj, vars.to.regress = c("nUMI", "percent.mito"), do.par=T, num.cores=10)

################################################################################

print('>>>CLUSTERING')

print('Finding variable genes...')
seurobj <- FindVariableGenes(seurobj, do.plot=F)
print('running pca')
seurobj <- RunPCA(seurobj, pcs.compute=50, do.print=F)

print('Running clustering...')
for (res in resolutions){
  seurobj <- FindClusters(seurobj, reduction.type = "pca", dims.use = 1:n.pcs, resolution = res, print.output = 0, save.SNN = TRUE)
}

################################################################################

print('>>>CALCULATING CELL CYCLE SCORES')

cc.genes <- readLines(cellcycle_genes)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
seurobj <- CellCycleScoring(seurobj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

################################################################################

print('>>>RUNNING TSNE')

seurobj <- RunTSNE(seurobj, reduction.use='pca', dims.use=1:n.pcs)

################################################################################

print('>>>SAVING DATASET')
print(paste('Saving dataset:', output))
saveRDS(seurobj, output)
