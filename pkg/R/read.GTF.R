read.GTF <- function(file, feature=c("exon", "CDS", "intron", "utr"), attributes=c("gene_id", "transcript_id", "gene_name", "transcript_name", "protein_id")){
  #require(GenomicRanges)
  cat("load GTF file ... \n")
  GTF <- read.table(file, sep="\t", header=FALSE, stringsAsFactors=FALSE, colClasses=c(rep("character", 3), rep("numeric", 2), rep("character", 4)))
  cat("parse attributes ... \n")
  attrs <- GTF[,9]
  attrs <- strsplit(attrs, split=" |;")
  for(i in 1:length(attrs)){
    ai <- match(attributes, attrs[[i]])
    if(sum(is.na(ai))==0){
      break
    }
  }
  Attrs <- do.call("cbind", lapply(attrs, function(x)x[ai+1]))
  Attrs <- t(Attrs)
  colnames(Attrs) <- attributes
  GR <- GRanges(seqnames=GTF[,1], IRanges(start=GTF[,4], end=GTF[,5]), strand=GTF[,7], source=GTF[,2], feature=GTF[,3], score=GTF[,6], frame=GTF[,8], Attrs)
  
  if(feature[1]=="all"){
    return(GR)
  }else{ #split by feature
    #exon
    if(length(grep("exon", feature, ignore.case=TRUE))>0){
      cat("extract exon ...\n")
      GRexon <- GR[GR@elementMetadata$feature=="exon"]
      GRL <- list(exon=sort(GRexon))
    }else{
      GRL <- list()
    }
    #CDS, stop_codon will be included
    if(length(grep("CDS", feature, ignore.case=TRUE))>0){
      cat("extract CDS ...\n")
      cds <- GR[GR@elementMetadata$feature %in% c("CDS", "stop_codon")]
      cds <- GRanges(seqnames=paste(seqnames(cds), cds@elementMetadata$transcript_id, sep="|"), ranges=ranges(cds), strand=strand(cds))
      GRcds <- reduce(cds)
      chrid <- do.call("rbind", strsplit(as.character(seqnames(GRcds)), split="\\|"))
      GRcds <- GRanges(seqnames=chrid[,1], ranges=ranges(GRcds), strand=strand(GRcds), transcript_id=chrid[,2])
      GRL <- c(GRL, CDS=GRcds)
    }
    if(length(grep("utr|intron", feature, ignore.case=TRUE))>0){
    #intron
    #define gaps
      cat("extract intron/utr ...\n")
      GRexon1 <- GRanges(seqnames=paste(seqnames(GRexon), GRexon@elementMetadata$transcript_id, sep="|"), ranges=ranges(GRexon), strand=strand(GRexon))
      GRgaps <- gaps(GRexon1)
      GRgaps <- GRgaps[start(GRgaps)!=1]
    #define intron, gaps
      chrid <- do.call("rbind", strsplit(as.character(seqnames(GRgaps)), split="\\|"))
      intron <- GRanges(seqnames=chrid[,1], ranges=ranges(GRgaps), strand=strand(GRgaps), transcript_id=chrid[,2])
      if(length(grep("intron", feature, ignore.case=TRUE))>0){
        GRL <- c(GRL, intron=intron)
      }
      if(length(grep("utr", feature, ignore.case=TRUE))>0){
        GRLcds <- split(GRcds, GRcds@elementMetadata$transcript_id)
        GRLexon <- split(GRexon, GRexon@elementMetadata$transcript_id)
        #filter
        GRLexon <- GRLexon[names(GRLexon) %in% names(GRLcds)]
        if(sum(names(GRLcds)!=names(GRLexon))>0)stop("exon and CDS not match")
        exons <- start(GRLexon)
        exone <- end(GRLexon)
        cdss <- start(GRLcds)
        cdse <- end(GRLcds)
        tids <- names(GRLexon)
        strands <- strand(GR)[match(tids, GR@elementMetadata$transcript_id)]
        seqs <- seqnames(GR)[match(tids, GR@elementMetadata$transcript_id)]
        
        u1 <- GRanges(seqs, IRanges(min(exons), min(cdss)-1), strand=strands, transcript_id=tids)
        u1 <- u1[start(u1)!=(end(u1)+1)]
        u2 <- GRanges(seqs, IRanges(max(cdse)+1, max(exone)), strand=strands, transcript_id=tids)
        u2 <- u2[(start(u2)-1)!=end(u2)]
        
        UTR1 <- diffGR(u1, intron)
        UTR2 <- diffGR(u2, intron)
        utr5 <- UTR1[strand(UTR1)=="+"]
        utr5 <- sort(c(utr5, UTR2[strand(UTR2)=="-"]))
        utr3 <- UTR1[strand(UTR1)=="-"]
        utr3 <- sort(c(utr3, UTR2[strand(UTR2)=="+"]))
        GRL <- c(GRL, utr3=utr3, utr5=utr5)
      }
      #GRL <- c(GRL, attributes=list(unique(Attrs)))
    }
    return(GRL)
  } 
}

write.GTF <- function(GR, file=""){
  GRtable <- cbind(seqname=as.character(seqnames(GR)), start=as.character(start(GR)), end=as.character(end(GR)), strand=as.character(strand(GR)), as.data.frame(GR@elementMetadata))
  idx <- match(c("seqname", "source", "feature", "start", "end", "score", "strand", "frame"), colnames(GRtable))
  idx1 <- na.omit(idx)
  attributes <- as.matrix(GRtable[,-idx1])
  atn <- colnames(attributes)

  atts <- c()
  for(i in 1:ncol(attributes)){
    atts1 <- paste(atn[i], attributes[,i])
    atts <- paste(atts, atts1, sep="; ")
  }
  atts <- sub("^;", "", atts)
  atts <- sub("\\s*\\w*\\sNA;?", "", atts)

  GRTable <- matrix(".", nrow=nrow(GRtable), ncol=9)
  if(is.null(na.action(idx1))){
    GRTable[,1:8] <- as.matrix(GRtable[,idx1])
  }else{
    GRTable[,c(1:8)[-na.action(idx1)]] <- as.matrix(GRtable[,idx1])
  }
  GRTable[,9] <- atts
  colnames(GRTable) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  write.table(GRTable, file=file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}
