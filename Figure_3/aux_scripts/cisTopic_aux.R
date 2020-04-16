# Get overlapping regions
getOverlapRegionsFromCoordinates <- function(
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