.libPaths('/home/cbmr/pytrik/libraries/')

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset.'),
  make_option(c('-o', '--outdir'), type='character', help='Output directory.'),
  make_option(c('-t', '--test'), type='character', help='Test to use.')
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

print('reading data')
data <- readRDS(opt$file)

output_prefix <- unlist(strsplit(opt$file, '/'))[length(unlist(strsplit(opt$file, '/')))]

#data <- SubsetData(data, cells.use=rownames(data@meta.data)[!data@meta.data$sample_name %in% 'Supra_4'])
#data <- SubsetData(data, cells.use=rownames(data@meta.data)[!data@meta.data$sample_name %in% 'Subq_4'])

print('setting type columns (brown or white)')
type <- unlist(lapply(data@meta.data$depot, function(x){
  if (x == 'Supra' || x == 'Peri'){
    return('brown')
  } else if (x == 'Subq' || x == 'Visce'){
    return('white')
  }
}))

#peri white
#supra white
#subq brown
#visce brown

sample2.white <- unlist(lapply(as.vector(data@meta.data$depot), function(x){
  if (x == 'Subq' || x == 'Visce'){
    return('white')
  } else {
    return(x)
  }
}))

sample2.brown <- unlist(lapply(as.vector(data@meta.data$depot), function(x){
  if (x == 'Supra' || x == 'Peri'){
    return('brown')
  } else {
    return(x)
  }
}))

data@meta.data['type'] <- type
data@meta.data['sample2.white'] <- sample2.white
data@meta.data['sample2.brown'] <- sample2.brown

#ident, group1, group2
steps <- list(c('sample2.white', 'Peri', 'white'),
              c('sample2.white', 'Supra', 'white'),
              c('sample2.brown', 'Visce', 'brown'),
              c('sample2.brown', 'Subq', 'brown'))

all.markers <- list()

for (x in steps){
  print(x)
  data <- SetAllIdent(data, id=x[1])
  print(2)
  markers <- c(1,2,3)
  print(3)
  print(x[2])
  print(x[3])
  markers <- FindMarkers(data, ident.1=x[2], ident.2=x[3], min.pct = 0.25, only.pos = F, thresh.use = 0.25, test.use = opt$test)
  print(4)
  all.markers[[paste(x[2], x[3], sep='.')]] = markers
}

#print(all.markers)

all.markers <- lapply(all.markers, function(x) cbind(x, gene=rownames(x)))
names(all.markers)
df.cluster_markers <- bind_rows(all.markers, .id='cluster')
df.cluster_markers <- df.cluster_markers[df.cluster_markers$p_val_adj < 0.05, ]

print('saving dataframe marker genes')
write.table(df.cluster_markers, file=paste(opt$outdir, '/', output_prefix, 'markergenes-crossed-depots', sep=''), sep='\t', row.names=F, quote=F)
