# Convert region name to data frame
regionName2DataFrame <- function(region_names){
  seqnames <- sapply(strsplit(region_names, split = ":"), "[", 1)
  coord <- sapply(strsplit(region_names, split = ":"), "[", 2)
  start <- sapply(strsplit(coord, split = "-"), "[", 1)
  end <- sapply(strsplit(coord, split = "-"), "[", 2)
  return(as.data.frame(cbind(seqnames, start, end)))
}

# Calculate distance index between a region and linked genes
distanceIndex2Gene <- function(region, TSS, data_bb){
  distance2genes <- distance(TSS, region)
  names(distance2genes) <- paste0(seqnames(TSS),':', start(TSS),'-', end(TSS))
  distance2genes <- distance2genes[-which(is.na(distance2genes))]
  distance2genes <- sort(distance2genes)
  region_name <- paste0(seqnames(region),':', start(region),'-', end(region))
  data_bb_region <- data_bb[which(as.vector(unlist(data_bb[,'sourceName'])) == region_name),]
  gene_index <- which(names(distance2genes) %in% unlist(as.vector(data_bb_region[,'targetName'])))
  return(gene_index)
}