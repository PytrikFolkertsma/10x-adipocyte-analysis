.libPaths('/home/cbmr/pytrik/libraries/')

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset.'),
  make_option(c('-o', '--outdir'), type='character', help='Output directory')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if (is.null(opt$file) || is.null(opt$outdir)){
  print_help(optionparser)
  quit()  
}

###################################################################################

library(Seurat)
library(parallel)
library(dplyr)

output_prefix <- unlist(strsplit(opt$file, '/'))[length(unlist(strsplit(opt$file, '/')))]
print('file name')
print(output_prefix)

print('Reading data')
data <- readRDS(opt$file)

samples <- as.character(data@meta.data$sample_name)
samples[which(data@meta.data$type == 'brown')] <- 'brown' 
data@meta.data['test'] <- samples
data <- SetAllIdent(data, id='test')

white_samples <- sort(unique(data@meta.data$sample_name[data@meta.data$depot %in% c('Visce', 'Subq')]))

print('Finding marker genes for white samples:')
print(white_samples)

cl <- makeCluster(8, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
list_of_dfs.all_markers <- parLapply(cl, white_samples, function(x) FindMarkers(data, ident.1 = x, ident.2='brown', min.pct = 0.1, test.use='negbinom'))
stopCluster(cl)

list_of_dfs.all_markers <- lapply(list_of_dfs.all_markers, function(x) cbind(x,gene=rownames(x))) # add gene name as column
names(list_of_dfs.all_markers) <- white_samples # set names so bind_rows will get the .id correct
df.cluster_markers <- bind_rows(list_of_dfs.all_markers, .id="sample_name") #combine list of dfs into a single data frame

print('saving dataframe marker genes')
write.table(df.cluster_markers, file=paste(opt$outdir, 'markers_', output_prefix, '_white-samples', sep=''), sep='\t', row.names=F, quote=F)
