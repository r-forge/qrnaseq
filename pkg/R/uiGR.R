#union or intersect by gene_id
uiGR <- function(GR, Reduce="union"){
  if(!("gene_id" %in% colnames(GR@elementMetadata))){
    stop("No gene_id found in elementMetadata")
  }
  GRL <- split(GR, GR@elementMetadata$gene_id)
  if(!is.na(pmatch(Reduce, "union"))){
    GRLu <- reduce(GRL)    
    GRu <- unlist(GRLu, use.names=FALSE)
    geneids <- rep(names(GRLu), lapply(GRLu, length))
    GRu@elementMetadata$gene_id <- geneids
    return(GRu)
  }else if(!is.na(pmatch(Reduce, "intersection"))){    
    GRa <- GRanges(seqnames=GR@elementMetadata$gene_id, ranges=ranges(GR), strand=strand(GR))
    cove <- coverage(GRa)
    covegr <- lapply(cove, function(x)tail(splitRanges(x),1)[[1]])

    coveg <- names(covegr)
    covei <- match(coveg, GR@elementMetadata$gene_id)
    coveseq <- seqnames(GR[covei])
    covestrand <- strand(GR[covei])
    coveL <- vector("list", length(cove))
    for(j in 1:length(cove)){
      coveL[[j]] <- GRanges(coveseq[j], covegr[[j]], covestrand[j], gene_id=rep(coveg[j], length(covegr[[j]])))
    }  
    GRi <- unlist(GRangesList(coveL))
    return(GRi)
  }else{
    stop("Wrong reduce method")
  }
}
