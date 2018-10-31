
### Runs CCA, aligns subspaces 30 CC's, runs tSNE on 15.
### Discarded cells are saved in data/alignment_discarded-cells

n.ccs <- 30 #Nr of CCs to compute 
n.ccs.use <- 15 #Nr of CCs to use for tSNE and clustering.

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the Seurat object to run the alignment on.'),
  make_option(c('-c', '--column'), type='character', help='Column to align on (e.g. sample_name or timepoint).')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if(is.null(opt$file) || is.null(opt$column)){
  print_help(optionparser)
  quit()
}

################################################################################

library(Rcpp)
library(magrittr)
library(Seurat)

################################################################################

print('LOADING DATA')

seurobj_all10x <- readRDS(opt$file)

################################################################################

print('>>>SPLITTING DATASET ON GIVEN COLUMN NAME')

all_groupnames <- unlist(seurobj_all10x@meta.data[opt$column])
groups <- unique(all_groupnames)
print('Subsets:')
print(groups)

objectlist <- list()
variable_genes <- c()

#for every group, create subset and find highly variable genes
for (i in c(1:length(groups))){
  print(paste(i, '- splitting data and finding variable genes for', groups[i]))
  objectlist[[i]] <- SubsetData(seurobj_all10x, cells.use=seurobj_all10x@cell.names[which(all_groupnames == groups[i])])
  objectlist[[i]] <- FindVariableGenes(objectlist[[i]], do.plot = FALSE)
  variable_genes <- c(variable_genes, rownames(head(objectlist[[i]]@hvg.info, n=1000)))
}

unique.variable_genes <- unique(variable_genes)
print(paste('Nr of variable genes:', length(unique.variable_genes)))

################################################################################

print('>>>RUNNING MULTICCA')

data <- RunMultiCCA(object.list = objectlist, genes.use=unique.variable_genes, num.ccs=n.ccs)

################################################################################

print('>>>ALIGNING SUBSPACES')

#Discard cells whose expression profile cannot be well-explained by low-dimensional CCA, 
#compared to low-dimensional PCA. Here cells are discarded when the ratio of the variance 
#explained by CCA is smaller than 0.5 (compared to the variance explained by PCA).

print('Discarding cells whose expression profile cannot be well explained by low-dim CCA compared to low-dim PCA...')
data <- CalcVarExpRatio(data, reduction.type = "pca", grouping.var = opt$column, dims.use = 1:10)
data.all.save <- data
data <- SubsetData(object = data, subset.name = "var.ratio.pca", accept.low = 0.5)
data.discard <- SubsetData(object = data.all.save, subset.name = "var.ratio.pca", accept.high = 0.5)

print(paste('Saving object discarded cells in output folder as ', opt$file, '-cca-discardedcells', sep=''))
saveRDS(data.discard, paste(opt$file, 'cca-discardedcells', sep='-'))

print(paste('Aligning subspaces on ', n.ccs, ' CCs...', sep=''))
data.aligned <- AlignSubspace(data, reduction.type = "cca", grouping.var = opt$column, dims.align = 1:n.ccs)

################################################################################

print('>>>RUNNING TSNE')
data.aligned <- RunTSNE(data.aligned, reduction.use='cca.aligned', dims.use=1:n.ccs.use)

print('>>>SAVING SEURAT OBJECT')
saveRDS(data.aligned, paste(opt$file, 'aligned', sep='-'))
