# .H <-  TFBSTools::shannon.entropy
.H <- function(pVect){
  # TODO: check negatives? 
  pVect <- pVect[pVect>0]
  - sum(pVect * log2(pVect))  # log or log2?
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}
calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
{
  if(class(AUC)=="aucellResults") AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}

# (points)
plotRSS <- function(rss, thisType, n=5)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,thisType], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    geom_point(color = "blue", size = 1) + 
    ggtitle(thisType) + 
    geom_label_repel(aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     na.rm=TRUE) +
    theme_classic()
}

plotRSS_heatmap <- function(rss, trh=NULL, row_names_gp=gpar(fontsize=5), order_rows=TRUE, cluster_rows=FALSE, name="RSS", ...)
{
  if(is.null(trh)) trh <- signif(quantile(rss, p=.97),2)
  
  library(ComplexHeatmap)
  rssSubset <- rss[rowSums(rss > trh)>0,]
  rssSubset <- rssSubset[,colSums(rssSubset > trh)>0]
  message("Showing regulons and cell types with any RSS > ", trh, " (dim: ", nrow(rssSubset), "x", ncol(rssSubset),")")
  
  if(order_rows)
  {
    maxVal <- apply(rssSubset, 1, which.max)
    rss_ordered <- rssSubset[0,]
    for(i in 1:ncol(rssSubset))
    {
      tmp <- rssSubset[which(maxVal==i),,drop=F]
      tmp <- tmp[order(tmp[,i], decreasing=FALSE),,drop=F]
      rss_ordered <- rbind(rss_ordered, tmp)
    }
    rssSubset <- rss_ordered
    cluster_rows=FALSE
  }
  
  Heatmap(rssSubset, name=name, row_names_gp=row_names_gp, cluster_rows=cluster_rows, ...)
} 
