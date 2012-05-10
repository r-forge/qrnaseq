#file="/user/songliu/u2/group/Qiang/projects/RNA-seq/data/annotation/knownGene.txt"
read.GenePred <- function(file, geneid=12, feature=c("exon", "CDS", "intron", "utr"), ...){
  gpred <- read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE, ...)
  nexon <- gpred[,8]
  exonStart <- strsplit(gpred[,9], split=",")
  exonEnd <- strsplit(gpred[,10], split=",")
  GRexon <- GRanges(seqnames=rep(gpred[,2], nexon), IRanges(as.numeric(unlist(exonStart))+1, as.numeric(unlist(exonEnd))+1), strand=rep(gpred[,3], nexon), transcript_id=rep(gpred[,1], nexon))
  if(!is.null(geneid)){
    GRexon@elementMetadata$gene_id <- rep(gpred[, geneid], nexon)
  }
  
  #exon
  if(length(grep("exon", feature, ignore.case=TRUE))>0){
    GRL <- list(exon=GRexon)
  }else{
    GRL <- list()
  }
  #define gaps
  GRexon1 <- GRanges(seqnames=paste(seqnames(GRexon), GRexon@elementMetadata$transcript_id, sep="|"), ranges=ranges(GRexon), strand=strand(GRexon))
  GRgaps <- gaps(GRexon1)
  GRgaps <- GRgaps[start(GRgaps)!=1]  
  #define intron
  chrid <- do.call("rbind", strsplit(as.character(seqnames(GRgaps)), split="\\|"))
  intron <- GRanges(seqnames=chrid[,1], ranges=ranges(GRgaps), strand=strand(GRgaps), transcript_id=chrid[,2])
  #define CDS
  if(length(grep("CDS", feature, ignore.case=TRUE))>0){
    gpred1 <- gpred[gpred[,6]!=gpred[,7],]
    gpredGRcds <- GRanges(gpred1[,2], IRanges(gpred1[,6]+1, gpred1[,7]), strand=gpred1[,3], transcript_id=gpred1[,1])
    GRcds <- diffGR(gpredGRcds, intron)
    GRL <- c(GRL, CDS=GRcds)
  }
  
  if(length(grep("intron", feature, ignore.case=TRUE))>0){
    GRL <- c(GRL, intron=intron)
  }

  #define utr
  if(length(grep("utr", feature, ignore.case=TRUE))>0){
    gpred1 <- gpred[gpred[,6]!=gpred[,7],]  
    u1 <- GRanges(gpred1[,2], IRanges(gpred1[,4]+1, gpred1[,6]), strand=gpred1[,3], transcript_id=gpred1[,1])
    u1 <- u1[start(u1)!=(end(u1)+1)]
    u2 <- GRanges(gpred1[,2], IRanges(gpred1[,7]+1, gpred1[,5]), strand=gpred1[,3], transcript_id=gpred1[,1])
    u2 <- u2[(start(u2)-1)!=end(u2)]
    UTR1 <- diffGR(u1, intron)
    UTR2 <- diffGR(u2, intron)
    utr5 <- UTR1[strand(UTR1)=="+"]
    utr5 <- sort(c(utr5, UTR2[strand(UTR2)=="-"]))
    utr3 <- UTR1[strand(UTR1)=="-"]
    utr3 <- sort(c(utr3, UTR2[strand(UTR2)=="+"]))
    GRL <- c(GRL, utr3=utr3, utr5=utr5)
  }
  return(GRL)
}
