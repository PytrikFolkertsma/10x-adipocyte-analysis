
library(Seurat, lib.loc='/home/cbmr/pytrik/libraries/')
library(parallel)
library(dplyr)

source('find-markers.R')

data <- readRDS('../output/10x-180831')
#data <- SubsetData(SetAllIdent(data, id='timepoint'), max.cells.per.ident = 1000)

T1T2T3 <- SubsetData(SetAllIdent(data, id='time_combined'), ident.use = 1)
T1T2T3 <- FindVariableGenes(T1T2T3, do.plot=F)
T1T2T3 <- RunPCA(T1T2T3, pcs.compute=10, do.print=F)
T1T2T3 <- FindClusters(T1T2T3, reduction.type = "pca", dims.use = 1:10, resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = T)

T4T5 <- SubsetData(SetAllIdent(data, id='time_combined'), ident.use = 2)
T4T5 <- FindVariableGenes(T4T5, do.plot=F)
T4T5 <- RunPCA(T4T5, pcs.compute=10, do.print=F)
T4T5 <- FindClusters(T4T5, reduction.type = "pca", dims.use = 1:10, resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = T)

markers_1 <- find_markers(T1T2T3, colname='res.1.5')
markers_2 <- find_markers(T4T5, colname='res.1.5')

write.table(markers_1, '../output/markergenes/markers_10x-180831-T1T2T3_res.1.5_negbinom', sep='\t', row.names=F, quote=F)
write.table(markers_2, '../output/markergenes/markers_10x-180831-T4T5_res.1.5_negbinom', sep='\t', row.names=F, quote=F)

genes <- unique(intersect(markers_1[markers_1$p_val_adj < 0.05, 'gene'], markers_2[markers_2$p_val_adj < 0.05, 'gene']))
print(paste('Nr of genes:', genes))

write.table(genes, file='../output/markergenes/10x-180831_genelist_T1T2T3-T4T5_res.1.5', sep='\n', row.names=F, quote=F)
