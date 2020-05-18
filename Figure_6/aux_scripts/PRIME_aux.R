readNheader <- function(matFile, motifAnnotFile="/media/data/lcb/icistarget/data/motifCollection/v9/motifs_count_md5_to_motif_names.tsv")
{
  ### Corresponding motif name:
  motifNames <- readLines(motifAnnotFile) # check if used v8 or v9
  motifNames <- motifNames[-1]
  motifNames <- sapply(motifNames, strsplit, "\t")
  table(lengths(motifNames))
  motifNames <- setNames(sapply(motifNames, function(x) x[2]), sapply(motifNames, function(x) x[1])) # choose the first motif name
  
  ## Read feather
  mat <- read_feather(matFile, columns = NULL)
  # Replace id by motif name
  colnames(mat)[2:ncol(mat)] <- unname(motifNames[colnames(mat)[2:ncol(mat)]])
  
  mat <- data.frame(mat, check.names=FALSE)
  row.names(mat) <- mat[,1]
  mat <- as.matrix(mat[,-1])
  
  if(any(is.na(colnames(mat))))  warning("missing motif names")
  return(mat)
}

# Differential enrichment with MAST
library(MAST)
DiffEnrichment <- function(deComparisons, Data){
  require(RcisTarget)
  data(motifAnnotations_dmel)
  mAnnot_data <- motifAnnotations_dmel[,1:2]
  mAnnot_data <- aggregate(mAnnot_data$TF, list(mAnnot_data$motif), paste, collapse=",")
  mAnnot <- as.vector(unlist(mAnnot_data[,2]))
  names(mAnnot) <- as.vector(unlist(mAnnot_data[,1]))
  Group <- Data[,1]
  names(Group) <- rownames(Data)
  Frame <- Data[,-1]
  for(DEcontrastName in colnames(deComparisons)) {
    DEclasses <- setNames(factor(deComparisons[, DEcontrastName]), rownames(deComparisons))
    flattenedEset <- melt(Frame)
    flattenedEset <- flattenedEset[, c(2,1,3)]
    colnames(flattenedEset) <- c("Motif", "Region", "exprs")
    flattenedEset <- data.frame(flattenedEset, DEclasses=DEclasses[as.character(flattenedEset$Region)])
    sca4DE <- FromFlatDF(dataframe=flattenedEset, 
                         idvars="Region", measurement="exprs", primerid="Motif", 
                         id=esetName, cellvars="DEclasses", class="SingleCellAssay")
    lrtOut <- LRT(sca4DE, "DEclasses", referent="rest")  # Each label vs "rest"
    lrtOut <- data.frame(lrtOut, pAdj=p.adjust(lrtOut$p.value, method="fdr"))
    lrtOut <- data.frame(lrtOut, sigLogP=-(log(lrtOut$pAdj) * lrtOut$direction))
    lrtOut <- lrtOut[order(lrtOut$sigLogP, decreasing=TRUE),]
    cellsSplit <- split(rownames(deComparisons), deComparisons[, DEcontrastName])
    regionsContrast <- which(rownames(Frame) %in% cellsSplit[[gsub("vsRest", "", DEcontrastName)]])
    rest <- which(rownames(Frame) %in% cellsSplit[['rest']])
    FC <- apply(logMat, 2, function(x) gtools::foldchange(mean(x[regionsContrast]),mean(x[rest])))
    lrtOut <- data.frame(lrtOut, signedFC=FC[as.character(lrtOut$primerid)])
    lrtOut <- lrtOut[order(lrtOut$sigLogP, decreasing=TRUE),]
    lrtOut <- cbind(lrtOut[,1:3], TF=as.vector(unlist(mAnnot[as.vector(unlist(lrtOut[[1]]$primerid))])), lrtOut[,4:ncol(lrtOut)])
    lrtOutTable[[DEcontrastName]] <- lrtOut
    colnames(lrtOutTable[[DEcontrastName]])[which(colnames(lrtOutTable[[DEcontrastName]]) %in% 'primerid')] <- 'motif'
    lrtOutTable[[DEcontrastName]] <- as.data.table(lrtOutTable[[DEcontrastName]])
    print(paste0(DEcontrastName, ' done!'))
    
  }
  return(lrtOutTable)
}


