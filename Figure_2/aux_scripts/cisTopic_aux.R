# library
library(GenomicRanges)
library(rtracklayer)
library(cisTopic)

# Select DR from cisTopicObject
.selectCoordinates <- function(
  object,
  target,
  method,
  dim=dim
){
  if (!target %in% c('cell', 'region')){
    stop('Please, provide target="cell" or "region".')
  } 
  
  if (method == 'tSNE'){
    if (is.null(object@dr[[target]][['tSNE']])){
      stop(paste0('Please, run first: cisTopicObject <- runtSNE(cisTopicObject, target="', target,'", ...).'))
    }
    
    coordinates <- object@dr[[target]][['tSNE']]
    
    if (ncol(coordinates) == 2 && dim > 2){
      stop(paste0('Please, run first: cisTopicObject <- runtSNE(cisTopicObject, target="', target,'", dim=3 ...)for a 3D tSNE.'))
    }
  }
  
  else if (method == 'Umap'){
    if (is.null(object@dr[[target]][['Umap']])){
      stop(paste0('Please, run first: cisTopicObject <- runUmap(cisTopicObject, target="', target,'", ...).'))
    }
    
    coordinates <- object@dr[[target]][['Umap']]
    
    if (ncol(coordinates) == 2 && dim > 2){
      stop(paste0('Please, run first: cisTopicObject <- runUmap(cisTopicObject, target="', target,'", n_components=3 ...)for a 3D tSNE.'))
    }
  }
  
  else if (method == 'PCA' || method == 'Biplot'){
    if (is.null(object@dr[[target]][['PCA']])){
      stop(paste0('Please, run first: cisTopicObject <- runPCA(cisTopicObject, target="', target,'", ...).'))
    }
    coordinates <- object@dr[[target]][['PCA']]$ind.coord
  }
  
  else if (method == 'DM'){
    if (is.null(object@dr[[target]][['DiffusionMap']])){
      stop(paste0('Please, run first: cisTopicObject <- runDM(cisTopicObject, target="', target,'", ...).'))
    }
    coordinates <- object@dr[[target]][['DiffusionMap']]
  }
  return(coordinates)
}

# Plot signatures or topic from cisTopicObject in red
SignatureEnrichmentBR <- function(object, target, method='Z-score', coordinates='tSNE', topic=NULL, signature=NULL, thrP=0.45, sort=TRUE){
  coordinates <- .selectCoordinates(object, target, coordinates, dim=2) 
  if (!is.null(topic)){
    modelMat <- modelMatSelection(object, target, method)
    modelMat <- t(apply(modelMat, 1, function(x) (x-min(x))/(max(x)-min(x))))
    modCols <- list(red=c(paste0('Topic', topic)))
    modCols <- lapply(modCols, function(x) sapply(x, function(topic) rownames(modelMat)[topic]))
    cellColChan <- sapply(modCols, function(modsCol) apply(as.matrix(modelMat[names(modsCol),]), 1, mean))
    cellColChan[is.na(cellColChan)] <- 0
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"],0,0, alpha=.8))
    names(cellCol) <- colnames(modelMat)
    if (sort == TRUE){
      coord_sorted <- names(sort(cellCol[rownames(coordinates)]))
      plot(coordinates[coord_sorted,], col=cellCol[coord_sorted], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE, main=paste('Topic', topic))
    }
    else{
      plot(coordinates, col=cellCol[rownames(coordinates)], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE, main=paste('Topic', topic))
    }
  }
  if (!is.null(signature)){
    if (target == 'cell'){
      modelMat <- t(cisTopicObject@cell.data[,signature,drop=F]) 
    } else {
      modelMat <- t(cisTopicObject@region.data[rownames(coordinates),signature,drop=F])
    }
    modelMat <- t(apply(modelMat, 1, function(x) (x-min(x))/(max(x)-min(x))))
    modCols <- list(red=signature)
    modCols <- lapply(modCols, function(x) sapply(x, function(topic) rownames(modelMat)[topic]))
    cellColChan <- sapply(modCols, function(modsCol) apply(as.matrix(modelMat[names(modsCol),]), 1, mean))
    cellColChan[is.na(cellColChan)] <- 0
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"],0,0, alpha=.8))
    names(cellCol) <- colnames(modelMat)
    if (sort == TRUE){
      coord_sorted <- names(sort(cellCol[rownames(coordinates)]))
      plot(coordinates[coord_sorted,], col=cellCol[coord_sorted], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE)
    }
    else{
      plot(coordinates, col=cellCol[rownames(coordinates)], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE)
    }
  }
}

# Retrieves the cisTopic object regions that intersect with the target regions (e.g. genes, etc...)
# Return is a list of region names (character) split by the @
# queryRegions:  GRanges
intersectToRegionSet <- function(cisTopicObject, targetRegions, splitBy=colnames(targetRegions@elementMetadata)[1], minOverlap=0.4)
{
  # To do: Check types
  if(!is.numeric(minOverlap) || (minOverlap<0 && minOverlap>=1)) stop("minOverlap should be a number between 0 and 1 (percentage of overlap between the regions).")
  if(is.null(targetRegions@elementMetadata) | ncol(targetRegions@elementMetadata)==0) targetRegions@elementMetadata <- DataFrame(name=as.character(targetRegions))
  if(!splitBy %in% colnames(targetRegions@elementMetadata)) stop("missing gene id annotation")
  if(length(table(lengths(targetRegions@elementMetadata[,splitBy]))) >1 ) stop("")
  
  # cisTopicObject
  cisTopicRegions <- cisTopicObject@region.ranges
  #elementMetadata(cisTopicRegions)[["regionNames"]] <- names(cisTopicRegions)
  
  cisTopicRegionsOverlap <- findOverlaps(cisTopicRegions, targetRegions,
                                         minoverlap=1,  #maxgap=0L, select="all",
                                         type="any", ignore.strand=TRUE)#, ...)
  
  if(minOverlap>0)
  {
    # In i-cisTarget, the default is 40% minimum overlap. Both ways: It takes the maximum percentage (of the peak or the ict region)
    # To reproduce those results:
    overlaps <- pintersect(cisTopicRegions[queryHits(cisTopicRegionsOverlap)], targetRegions[subjectHits(cisTopicRegionsOverlap)])
    percentOverlapHuman <- width(overlaps) / width(cisTopicRegions[queryHits(cisTopicRegionsOverlap)])
    percentOverlapPeaks <- width(overlaps) / width(targetRegions[subjectHits(cisTopicRegionsOverlap)])
    maxOverlap <- apply(cbind(percentOverlapHuman, percentOverlapPeaks), 1, max)
    cisTopicRegionsOverlap <- cisTopicRegionsOverlap[maxOverlap > minOverlap]
  }
  
  hitsPerGene <- split(unname(as.character(cisTopicRegions[queryHits(cisTopicRegionsOverlap)])), unlist(targetRegions[subjectHits(cisTopicRegionsOverlap)]@elementMetadata[,splitBy]))
  regionsPerGene <- lapply(hitsPerGene, unique)
  
  missingGenes <- targetRegions@elementMetadata[,splitBy][which(!targetRegions@elementMetadata[,splitBy] %in% names(regionsPerGene))]
  regionsPerGene <- c(missing=list(missingGenes), regionsPerGene)
  
  return(regionsPerGene)
}

# Overlap
.getSignatureRegions <- function(
  object,
  signature,
  minOverlap=0.4,
  ...
){
  if (!is.data.frame(signature)){
    regions <- read.table(signature)
    colnames(regions)[1:3] <- c('seqnames', 'start', 'end')
    regions <- makeGRangesFromDataFrame(regions)
  }
  else if (is.data.frame(signature)){
    regions <- makeGRangesFromDataFrame(signature)
  }
  else{
    stop('The signature is not an existing file or a dataframe.')
  }
  
  coordinates <- object@region.ranges
  regionsSignature <- .getOverlapRegionsFromCoordinates(coordinates, regions, minOverlap=minOverlap, ...)
  message(paste('The signature contains', length(regions), 'of which', nrow(as.data.frame(regionsSignature)), 'overlap the regions in the set.'))
  return(regionsSignature)
}

.getOverlapRegionsFromCoordinates <- function(
  coordinates,
  regions,
  minOverlap=0.4,
  overlapping=TRUE,
  ...)
{
  dbRegionsOverlap <- findOverlaps(regions, coordinates, type='any', select="all", ignore.strand=TRUE, ...)
  
  if(minOverlap>0){
    overlaps <- pintersect(regions[queryHits(dbRegionsOverlap)], coordinates[subjectHits(dbRegionsOverlap)])
    percentOverlapCoordinates <- width(overlaps) / width(regions[queryHits(dbRegionsOverlap)])
    percentOverlapRegions <- width(overlaps) / width(coordinates[subjectHits(dbRegionsOverlap)])
    maxOverlap <- apply(cbind(percentOverlapCoordinates, percentOverlapRegions), 1, max)
    dbRegionsOverlap <- dbRegionsOverlap[maxOverlap > minOverlap]
    maxOverlap <- maxOverlap[which(maxOverlap > minOverlap)]
  }
  
  selectedRegions <- regions[queryHits(dbRegionsOverlap)]
  selectedRegions <- paste(as.vector(seqnames(selectedRegions)), ':', as.vector(start(selectedRegions)), '-', as.vector(end(selectedRegions)), sep='')
  selectedCoordinates <- names(coordinates[subjectHits(dbRegionsOverlap)])
  
  selectedMapped <- data.frame(selectedCoordinates, selectedRegions, maxOverlap, row.names=NULL)
  
  if (overlapping != TRUE){
    if(any(duplicated(selectedRegions))){
      selectedMapped <- selectedMapped[order(as.vector(selectedMapped$selectedRegions), -abs(selectedMapped$maxOverlap)), ]
      selectedMapped <- selectedMapped[!duplicated(as.vector(selectedMapped$selectedRegions)), ]
    }
  }
  return(selectedMapped)
}

# Get -log10(qval) from MACS output into cisTopicObject

addMACSlogqval <- function(object, signature, label){
  selectedMapped <- .getSignatureRegions(object, signature, overlapping=FALSE)
  selectedMapped <- selectedMapped[order(as.vector(selectedMapped$selectedCoordinates), -abs(selectedMapped$maxOverlap)), ]
  regionsSignatures <- selectedMapped[!duplicated(as.vector(selectedMapped$selectedCoordinates)), ]
  data_diff <- read.table(signature)
  logqval <- data_diff[,9,drop=FALSE] 
  rownames(logqval) <- paste0(data_diff[,1], ':', data_diff[,2], '-', data_diff[,3])
  regionsSignatures <- cbind(regionsSignatures, logqval[regionsSignatures[,1],])
  realSignature <- regionsSignatures[,c(1,4)]
  rownames(realSignature) <- realSignature[,1]
  oneqval <- rep(min(realSignature[,2]), nrow(object@region.data))
  names(oneqval) <- rownames(object@region.data)
  oneqval[rownames(realSignature)] <- realSignature[,2]
  oneqval <- as.data.frame(log(oneqval))
  colnames(oneqval) <- paste0('-log(qval)_', label)
  rownames(oneqval) <- rownames(object@region.data)
  object <- addRegionMetadata(object, oneqval)
  return(object)
}