
data <- readRDS('output/10x-180831')

metadata <- data@meta.data
metadata['rownames'] <- rownames(metadata)

#add dummy variables to metadata for timepoint+depot combinations
depots_timepoints <- unlist(apply(data@meta.data[,c('depot', 'timepoint')], 1, function(x){
  if (x[['timepoint']] %in% c('T1', 'T2', 'T3')){
    return(paste(x[['depot']], 'T1T2T3', sep='.'))
  } else {
    return(paste(x[['depot']], x[['timepoint']], sep='.'))
  }
}))
metadata['depots_timepoints'] <- depots_timepoints
metadata <- metadata %>% mutate(value=1) %>% spread(depots_timepoints, value, fill=0)

#add dummy variables for fat type + time
type_timepoints <- unlist(apply(data@meta.data[,c('type', 'timepoint')], 1, function(x){
  #if (x[['timepoint']] %in% c('T1', 'T2', 'T3')){
  #  return(paste(x[['type']], 'T1T2T3', sep='.'))
  #} else {
    return(paste(x[['type']], x[['timepoint']], sep='.'))
  #}
}))
metadata['type_timepoints'] <- type_timepoints
metadata <- metadata %>% mutate(value=1) %>% spread(type_timepoints, value, fill=0)

#replace NA's
metadata[which(is.na(data@meta.data$State)),'State'] <- 1

#add dummy variables for Monocle state
states <- unlist(lapply(metadata$State, function(x){
  if(x == 1){
    return('pre_branch')
  } else if (x==2){
    return('brown_branch')
  } else if (x==3) {
    return('white_branch')
  }
}))

metadata['states'] <- states
metadata <- metadata %>% mutate(value=1) %>% spread(states, value, fill=0)

state_timepoints <- unlist(apply(data@meta.data[,c('State', 'timepoint')], 1, function(x){
  if (x[['State']] == 1){
    return(paste('pre_branch', x[['timepoint']], sep='.'))
  } else if (x[['State']] == 2){
    return(paste('brown_branch', x[['timepoint']], sep='.'))
  } else {
    return(paste('white_branch', x[['timepoint']], sep='.'))
  }
}))

metadata['state_timepoints'] <- state_timepoints
metadata <- metadata %>% mutate(value=1) %>% spread(state_timepoints, value, fill=0)

#add dummy variables for timepoint -> are any modules related to general development?
metadata['timepoints'] <- unlist(lapply(as.vector(data@meta.data$timepoint), function(x){
  return(as.integer(substring(x, 2, 2)))
}))


#add columns
metadata['depots_timepoints'] <- depots_timepoints
metadata['type_timepoints'] <- type_timepoints
metadata['states'] <- states


#metadata_corr_col
columns <- colnames(metadata)
columns <- columns[match('Peri.T1T2T3', columns):match('timepoints', columns)]
dput(columns)

rownames(metadata) <- metadata$rownames
metadata$rownames <- NULL
data <- AddMetaData(data, metadata[columns])

#c("Peri.T1T2T3", "Peri.T4", "Peri.T5", "Subq.T1T2T3", "Subq.T4",
"Subq.T5", "Supra.T1T2T3", "Supra.T4", "Supra.T5", "Visce.T1T2T3",
"Visce.T4", "Visce.T5", "brown.T1", "brown.T2", "brown.T3", "brown.T4",
"brown.T5", "white.T1", "white.T2", "white.T3", "white.T4", "white.T5",
"brown_branch", "pre_branch", "white_branch", "brown_branch.T4",
"brown_branch.T5", "pre_branch.T1", "pre_branch.T2", "pre_branch.T3",
"pre_branch.T4", "pre_branch.T5", "white_branch.T2", "white_branch.T3",
"white_branch.T4", "white_branch.T5", "timepoints")
save(data, file='output/10x-180831.RData')

