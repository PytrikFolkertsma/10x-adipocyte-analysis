library(gprofiler)

df.cluster_markers <- read.table('output/markergenes/180504/markers_10x-180504_res.0.5_negbinom', header=T, sep='\t')
df.cluster_markers <- df.cluster_markers[df.cluster_markers$p_val_adj < 0.05,]
df.cluster_markers <- df.cluster_markers[df.cluster_markers$cluster == 12,]
df.cluster_markers <- df.cluster_markers[order(-df.cluster_markers$avg_logFC),]
mixture_pos_markers <- df.cluster_markers[df.cluster_markers$avg_logFC > 0,]
mixture_neg_markers <- df.cluster_markers[df.cluster_markers$avg_logFC < 0,]

gsea <- gprofiler(query=as.vector(mixture_pos_markers$gene[1:51]), organism='hsapiens', significant=T, ordered_query = F, src_filter='GO')
gsea <- gsea[order(gsea$p.value),]

write.table(gsea, file='tables/GSEA_mixture-cluster.txt', sep='\t', quote=F)
write.table(data.frame(gsea$term.id, gsea$p.value), file='tables/Revigo_input_mixture-cluster.txt', sep='\t', quote=F, row.names=F)

##########################################################################################

revigo <- read.table('tables/Revigo_mixture_cluster.csv', sep=',', header=T)
revigo <- revigo[revigo$eliminated == 0,]

revigo_merged <- merge(x=data.frame(revigo$term_ID, revigo$description, revigo$log10.p.value, revigo$frequency), 
                       y=data.frame(gsea$query.size, gsea$term.id, gsea$term.size, gsea$overlap.size, gsea$p.value, gsea$intersection), 
                       by.x='revigo.term_ID', 
                       by.y='gsea.term.id')

revigo_merged <- revigo_merged[order(revigo_merged$gsea.p.value),]

write.table(revigo_merged, 'tables/Revigo_mixture-cluster_cleaned.txt', sep='\t', quote=T, row.names=F)
