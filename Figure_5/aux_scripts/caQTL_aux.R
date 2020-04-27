# Calculate snps with a delta score above a threshold

sumAbove <- function(vector, thr=3){ 
  vector[abs(vector) <= thr] <- 0
  vector[abs(vector) > thr] <- 1
  return(sum(vector, na.rm=TRUE))
}

# Calculate aggregated delta score above a threshold

sumAbove_aggr <- function(vector, thr=3){
  vector[abs(vector) <= thr] <- 0
  return(sum(vector, na.rm=TRUE))
}

# Heatmap

dotheatmap <- function (enrichmentDf,
                        # top = 5, order.by = "Binom_Fold_Enrichment",
                        var.x="Topic", var.y="ID", 
                        var.col="FC", col.low="brown1", col.mid="floralwhite", col.high="dodgerblue", 
                        var.size="p.adjust", min.size=1, max.size=8
)
{
  require(data.table)
  require(ggplot2)
  
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  limit <- max(abs(enrichmentDf[,var.col])) * c(-1, 1)
  p <- ggplot(data=enrichmentDf, mapping=aes_string(x=var.x, y=var.y)) + 
    geom_point(mapping=aes_string(size=var.size, color=var.col)) +
    scale_radius(range=c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10), limit = limit) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=90, hjust=1))
  print(p)
}

# SNP name to genomic ranges 

SNPNames2GRanges <- function(SNP_names){
  data <- list()
  for (SNP_name in SNP_names){
    seqnames <- sapply(strsplit(SNP_name, split = ":"), "[", 1)
    coord <- sapply(strsplit(SNP_name, split = ":"), "[", 2)
    start <- as.numeric(as.vector(unlist(sapply(strsplit(coord, split = "__"), "[", 1))))
    end <- start+1
    data[[SNP_name]] <- as.data.frame(cbind(seqnames, start, end))
  }
  data <- rbindlist(data)
  return(makeGRangesFromDataFrame(data))
}

# Color by region accessibility

# Plot signatures or topic from cisTopicObject in red
RegionEnrichmentGreen <- function(object, pred.dist, region, target, method='Z-score', coordinates='tSNE', thrP=1, sort=TRUE){
  coordinates <- .selectCoordinates(object, target, coordinates, dim=2) 
  modelMat <- pred.dist[region,,drop=F]
  modelMat <- t(apply(modelMat, 1, function(x) (x-min(x))/(max(x)-min(x))))
  modCols <- list(red=region)
  modCols <- lapply(modCols, function(x) sapply(x, function(topic) rownames(modelMat)[topic]))
  cellColChan <- sapply(modCols, function(modsCol) apply(as.matrix(modelMat[names(modsCol),]), 1, mean))
  cellColChan[is.na(cellColChan)] <- 0
  cellColChan[cellColChan > thrP] <- thrP
  cellColChan <- (cellColChan-min(cellColChan))/(max(cellColChan)-min(cellColChan))
  cellCol <- apply(cellColChan, 1, function(x) rgb(0,x["red"],0, alpha=.8))
  names(cellCol) <- cisTopicObject@cell.names
  if (sort == TRUE){
    coord_sorted <- names(sort(cellCol[rownames(coordinates)]))
    plot(coordinates[coord_sorted,], col=cellCol[coord_sorted], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE)
  }
  else{
    plot(coordinates, col=cellCol[rownames(coordinates)], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE)
  }
}

# Select coordinates

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
