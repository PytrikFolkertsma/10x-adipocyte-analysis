library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset.'),
  make_option(c('-c', '--colname'), type='character', help='Metadata column to use for finding marker genes. Will test each group against the rest.'),
  make_option(c('-o', '--outdir'), type='character', help='Output directory.'),
  make_option(c('-t', '--test'), type='character', help='Test to use. Default is negative binomial.', default='negbinom')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if (is.null(opt$file) || is.null(opt$colname) || is.null(opt$outdir)){
  print_help(optionparser)
  quit()  
}

###################################################################################

library(Seurat)
library(parallel)
library(dplyr)

output_prefix <- unlist(strsplit(opt$file, '/'))[length(unlist(strsplit(opt$file, '/')))]
print(paste('file name:', output_prefix))

print('Reading data')
data <- readRDS(opt$file)
data <- SetAllIdent(data, id=opt$colname)

print('Finding marker genes')
#Identify cluster markers
cluster_ids <- sort(unique(data@meta.data[,opt$colname])) # numeric vector, 0, 1, 2, 3, ... (N_CLUSTERS-1)
print('cluster ids:')
print(cluster_ids)

cl <- makeCluster(8, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
list_of_dfs.all_markers <- parLapply(cl, cluster_ids, function(x) FindMarkers(data, ident.1 = x, min.pct = 0.1, test.use=opt$test))
stopCluster(cl)

list_of_dfs.all_markers <- lapply(list_of_dfs.all_markers, function(x) cbind(x,gene=rownames(x))) # add gene name as column
names(list_of_dfs.all_markers) <- cluster_ids # set names so bind_rows will get the .id correct
df.cluster_markers <- bind_rows(list_of_dfs.all_markers, .id="cluster") #combine list of dfs into a single data frame

print('saving dataframe marker genes')
write.table(df.cluster_markers, file=paste(opt$outdir, 'markers_', output_prefix, '_', opt$colname, '_', opt$test, sep=''), sep='\t', row.names=F, quote=F)
