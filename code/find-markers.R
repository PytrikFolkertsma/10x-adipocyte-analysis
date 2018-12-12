

find_markers <- function(data, colname, test='negbinom', num_nodes=8){
  print('Finding marker genes')
  cluster_ids <- sort(unique(data@meta.data[,colname]))
  print(paste('cluster ids:', cluster_ids))

  cl <- makeCluster(num_nodes, type = "FORK")
  clusterEvalQ(cl, library(Seurat, lib.loc='/home/cbmr/pytrik/libraries/'))
 #clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(rlang))
  list_of_dfs.all_markers <- parLapply(cl, cluster_ids, function(x) FindMarkers(data, ident.1 = x, min.pct = 0.1, test.use=test))
  stopCluster(cl)

  list_of_dfs.all_markers <- lapply(list_of_dfs.all_markers, function(x) cbind(x,gene=rownames(x)))
  names(list_of_dfs.all_markers) <- cluster_ids
  df.cluster_markers <- bind_rows(list_of_dfs.all_markers, .id="cluster")
  print('NAMES')
  print(names(df.cluster_markers))
  print(head(df.cluster_markers))
  df.cluster_markers <- df.cluster_markers[df.cluster_markers$p_val_adj < 0.05, ]
  return(df.cluster_markers)
}

if (sys.nframe() == 0){

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

  ########################################################################

  library(Seurat, lib.loc = '/home/cbmr/pytrik/libraries/')
  library(parallel)
  library(dplyr)

  output_prefix <- unlist(strsplit(opt$file, '/'))[length(unlist(strsplit(opt$file, '/')))]
  print(paste('file name:', output_prefix))

  print('Reading data')
  data <- readRDS(opt$file)
  print('test')
  data <- SetAllIdent(data, id=opt$colname)
  markers <- find_markers(data, colname=opt$colname, test=opt$test)

  print('saving dataframe marker genes')
  write.table(markers, file=paste(opt$outdir, 'markers_', output_prefix, '_', opt$colname, '_', opt$test, sep=''), sep='\t', row.names=F, quote=F)
}
