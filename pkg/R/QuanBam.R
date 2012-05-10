QuanBam <- function(bam, GR, method="both", Reduce="union", ovtype="within", revstrand=FALSE, ingene=FALSE, ...){
  #require(ShortRead)
  Rreads <- readBamGappedAlignments(bam, use.names=TRUE, ...)
  cat(length(Rreads), "reads are loaded\n")
#reduce
  if(!is.null(Reduce)){
    GR <- uiGR(GR, Reduce=Reduce)
    cat(paste(length(GR), "ranges left after reduced\n"))
  }
  
#filter1: match with first segment, type: within
  wids <- sub("M.*", "", cigar(Rreads))
  Rreadsf <- GRanges(seqnames=seqnames(Rreads), IRanges(start=start(Rreads), width=as.numeric(wids)))
  ov1 <- suppressWarnings(findOverlaps(Rreadsf, GR, type=ovtype))
  if("matchMatrix" %in% slotNames(ov1)){
    ov1m <- ov1@matchMatrix
  }else if("queryHits" %in% slotNames(ov1)){
    ov1m <- cbind(ov1@queryHits, ov1@subjectHits)
  }

  strandf <- as.character(strand(Rreads))
  nmf <- names(Rreads)

  cat("Filtering:\n")
  filter0 <- unique(ov1m[,1])
  filter1 <- rep(FALSE, length(nmf))
  filter1[filter0] <- TRUE
  cat(paste(sum(!filter1), " reads are not '", ovtype, "' mapped\n", sep=""))
  
#filter2: opposite direction  
  if(revstrand){
    nsf <- paste(nmf, strandf)
    dnsf <- nsf[duplicated(nsf)]
    samd <- sub(" .", "", dnsf )
    filter2 <- !(nmf %in% samd)
    cat(paste(sum(!filter2), "reads are not in reverse direction between pairs\n"))
  }else{
    filter2 <- TRUE
  }

#filter3: not in the same gene
  if(ingene){
    nmf3 <- nmf[ov1m[,1]]
    gids3 <- GR@elementMetadata$gene_id[ov1m[,2]]
    dnmf3 <- unique(nmf3[duplicated(nmf3)])
    n2g <- paste(nmf3, gids3)[nmf3 %in% dnmf3]
    dn2g <- n2g[duplicated(n2g)]
    samg <- sub(" ENSG.*", "", dn2g)
    filter3 <- !(nmf %in% setdiff(dnmf3, samg))
    cat(paste(sum(!filter3), "reads are not in the same genes between pairs\n"))
  }else{
    filter3 <- TRUE
  }

  RreadsF <- Rreads[filter1 & filter2 & filter3]
  Rreadsf1 <- Rreadsf[filter1 & filter2 & filter3]
  cat(paste(length(RreadsF), "reads are left\n"))

#Count
  if(method=="count" | method=="both"){
    cat("Summarize counts ...\n")
    counts <- suppressWarnings(countOverlaps(GR, Rreadsf1, type="any"))
    GR@elementMetadata$counts <- counts
    #RPKM
  }

#depth
  if(method=="maxdepth" | method=="both"){
    cat("Summarize maximum depths ...\n")
    cove <- coverage(RreadsF)
    GRL <- split(GR, seqnames(GR))
    GRL <- GRL[unlist(lapply(GRL, function(x)length(x)>0))]
    depthl <- vector(mode="list", length(GR))
    for(i in 1:length(GRL)){
      chi <- match(names(GRL)[i], names(cove))
      if(is.na(chi)){
        depthl[[i]] <- rep(0, length(GRL[[i]]))
      }else{
        r1 <- cbind(start(GRL[[i]]), end(GRL[[i]]))
        depthl[[i]] <- apply(r1, 1, function(x)max(cove[[chi]][x[1]:x[2]]))
      }
    }
    posl <- lapply(GRL, function(x)paste(seqnames(x), start(x), end(x)))
    posL <- unlist(posl)
    pos <- paste(seqnames(GR), start(GR), end(GR))
    posi <- match(pos, posL)
    depths <- unlist(depthl)[posi]
    GR@elementMetadata$mdepths <- depths
  }
  return(GR)
}
