
visce_T1T2T3 <- SubsetData(SetAllIdent(visce, id='time_combined'), ident.use=1)
visce_T4T5 <- SubsetData(SetAllIdent(visce, id='time_combined'), ident.use=2)

visce_T1T2T3 <- FindVariableGenes(visce_T1T2T3, do.plot=F)
visce_T1T2T3 <- RunPCA(visce_T1T2T3)
#PCElbowPlot(visce_T1T2T3)
visce_T1T2T3 <- FindClusters(visce_T1T2T3, reduction.type='pca', dims.use=1:10, resolution=1, force.recalc = T)
visce_T1T2T3 <- RunTSNE(visce_T1T2T3, reduction.use='pca', dims.use=1:10)
TSNEPlot(visce_T1T2T3, group.by='res.1', pt.size=0.1)
markers_visce_T1T2T3 <- FindAllMarkers(SetAllIdent(visce_T1T2T3, id='res.1'), test.use='negbinom')
markers_visce_T1T2T3_filtered <- markers_visce_T1T2T3[markers_visce_T1T2T3$p_val_adj < 0.05,]

visce_T4T5 <- FindVariableGenes(visce_T4T5, do.plot=F)
visce_T4T5 <- RunPCA(visce_T4T5)
#PCElbowPlot(visce_T1T2T3)
visce_T4T5 <- FindClusters(visce_T4T5, reduction.type='pca', dims.use=1:10, resolution=1, force.recalc = T)
visce_T4T5 <- RunTSNE(visce_T4T5, reduction.use='pca', dims.use=1:10)
TSNEPlot(visce_T4T5, group.by='res.1', pt.size=0.1)
markers_visce_T4T5 <- FindAllMarkers(SetAllIdent(visce_T4T5, id='res.1'), test.use='negbinom')
markers_visce_T4T5_filtered <- markers_visce_T4T5[markers_visce_T4T5$p_val_adj < 0.05,]

genes_visce <- unique(c(markers_visce_T1T2T3_filtered$gene, markers_visce_T4T5_filtered$gene))

write.table(genes_visce, col.names=F, row.names=F, file='output/markergenes/180831/monocle_genelist_visce_T1T2T3_T4T5_res.1', sep='\n', quote=F)



write.table(unique(markers_visce_T4T5_filtered$gene), col.names=F, row.names=F, file='output/markergenes/180831/monocle_genelist_visce_T4T5_res.1', sep='\n', quote=F)
write.table(unique(markers_subq_T4T5_filtered$gene), col.names=F, row.names=F, file='output/markergenes/180831/monocle_genelist_subq_T4T5_res.1', sep='\n', quote=F)
write.table(unique(markers_peri_T4T5_filtered$gene), col.names=F, row.names=F, file='output/markergenes/180831/monocle_genelist_peri_T4T5_res.1', sep='\n', quote=F)
write.table(unique(markers_supra_T4T5_filtered$gene), col.names=F, row.names=F, file='output/markergenes/180831/monocle_genelist_supra_T4T5_res.1', sep='\n', quote=F)





