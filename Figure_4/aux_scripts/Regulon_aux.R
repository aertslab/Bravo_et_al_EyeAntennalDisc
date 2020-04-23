# Read gmt
read_gmt <- function(file.gmt){
  x0 <- scan(file.gmt, what="", sep="\n", quiet=TRUE)
  x1 <- strsplit(x0, "\t")
  x2 <- lapply(x1, FUN=function(v){
    list(name=v[1], description=v[2], genes=v[c(-1,-2)])
  })
  names(x2) <- lapply(x2, FUN=function(vv) vv$name)
  x3 <- lapply(x2, FUN=function(vv) vv$name)
  return(x2)
}

# Extract region linked to TF motif from cisTarget results

linkRegulons <- function(motifEnrichment, data_bb, modules){
  # Rename motifDb
  names(motifEnrichment) <- paste0(1:length(names(motifEnrichment)), '__', names(motifEnrichment))
  names(modules) <-  paste0(1:length(names(modules)), '__', names(modules))
  motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
    cbind(motifDb=dbName, motifEnrichment[[dbName]])
  }))
  motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
  motifEnrichment_selfMotifs$highlightedTFs <- motifEnrichment_selfMotifs$motifDb
  # Reformat
  motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs, 1, function(oneMotifRow) {
    genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
    oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
    data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
  })
  motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
  colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
  motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)
  # Reformat
  regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
    # print(unique(tfTargets$TF))
    tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
      highConfAnnot <- "**" %in% enrOneGene$annot
      enrOneGeneByAnnot <- enrOneGene
      if(highConfAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
      bestMotif <- which.max(enrOneGeneByAnnot$NES)
      cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene),
            bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]),
            highConfAnnot=highConfAnnot)
    })), stringsAsFactors=FALSE)
    tfTable[order(tfTable$NES, decreasing = TRUE),]
  })
  regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
  colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "highConfAnnot")
  # Create regulons
  regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$highConfAnnot)
  regulons <- NULL
  if(!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]]))
  {
    regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
  }
  regulons_extended <- NULL
  if(!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]]))
  {
    regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(unlist(x[,"gene"])))
    regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), names(regulons_extended))
    names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
  }
  regulons <- c(regulons, regulons_extended)
  regulons <- sapply(regulons, function(x) gsub('dmel_r6.02__', '', x))
  # Regions to genes
  gene_regulons <- lapply(regulons, function(x) data_bb[which(as.vector(unlist(data_bb[,'sourceName'])) %in% x), ])
  gene_regulons <- lapply(1:length(gene_regulons), function(i) gene_regulons[[i]][which(as.vector(unlist(gene_regulons[[i]][,'name'])) %in% modules[[gsub('_extended', '',names(gene_regulons)[i])]]), ])
  names(gene_regulons) <- sapply(strsplit(names(regulons), split = "__"), "[", 2)
  gene_regulons_merge <- lapply(unique(names(gene_regulons)), function(x) rbindlist(gene_regulons[which(names(gene_regulons) %in% x)]))
  gene_regulons_merge <- lapply(gene_regulons_merge, function(x) x[!duplicated(x)])
  names(gene_regulons_merge) <- unique(names(gene_regulons))
  return(gene_regulons_merge)
}

# Plot enrichment


plotEnrichmentCustom <- function(pathway, stats, gseaParam = 1, v=120){
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = "dodgerblue", 
                                                      size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
                                                                               linetype = "dashed") +  geom_vline(xintercept = v, colour = "grey", 
                                                                                                                  linetype = "dashed") + geom_hline(yintercept = min(bottoms), 
                                                                                                                                                    colour = "red", linetype = "dashed") + geom_hline(yintercept = 0, 
                                                                                                                                                                                                      colour = "black") + geom_line(color = "dodgerblue") + theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
                                                               y = -diff/2, xend = x, yend = diff/2), size = 0.2) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "Rank", y = "Enrichment Score")
  g
}
