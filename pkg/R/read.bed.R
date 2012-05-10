#file="/user/songliu/u2/group/Qiang/projects/RNA-seq/data/annotation/Homo_sapiens.GRCh37.65.clean.bed"

read.bed <- function(file, feature=c("exon", "CDS", "intron", "utr")){
  bed <- read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  if(ncol(bed)<12){stop("Additional fields are required")}
  colnames(bed) <- c("chrom", "chromStart", "chromEnd", "transcript_id", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
  bsizes <- strsplit(bed[,11], split=",")
  bstarts <- strsplit(bed[,12], split=",")
  bedGRexon <- GRanges(bed[,1], IRanges(bed[,2]+1, bed[,3]), strand=bed[,6], bed[, -c(1,2,3,6,11,12)])
  bedGR1 <- rep(bedGRexon, bed[,10])
  GRexon <- GRanges(seqnames(bedGR1), IRanges(start=(start(bedGR1)+as.numeric(unlist(bstarts))), width=as.numeric(unlist(bsizes))), strand(bedGR1), bedGR1@elementMetadata)
  if(length(grep("exon", feature, ignore.case=TRUE))>0){
    GRL <- list(exon=sort(GRexon))
  }else{
    GRL <- list()
  }
  #define gaps
  GRexon1 <- GRanges(seqnames=paste(seqnames(GRexon), GRexon@elementMetadata$transcript_id, sep="|"), ranges=ranges(GRexon), strand=strand(GRexon))
  GRgaps <- gaps(GRexon1)
  GRgaps <- GRgaps[start(GRgaps)!=1]
  #define intron, gaps
  chrid <- do.call("rbind", strsplit(as.character(seqnames(GRgaps)), split="\\|"))
  intron <- GRanges(seqnames=chrid[,1], ranges=ranges(GRgaps), strand=strand(GRgaps), transcript_id=chrid[,2])

  #define CDS
  if(length(grep("CDS", feature, ignore.case=TRUE))>0){
    bed1 <- bed[bed[,7]!=bed[,8],]
    bedGRcds <- GRanges(bed1[,1], IRanges(bed1[,7]+1, bed1[,8]), strand=bed1[,6], transcript_id=bed1[,4])
    GRcds <- diffGR(bedGRcds, intron)
    GRL <- c(GRL, CDS=sort(GRcds))
  }
  
  if(length(grep("intron", feature, ignore.case=TRUE))>0){
    GRL <- c(GRL, intron=intron)
  }
  #define utr
  if(length(grep("utr", feature, ignore.case=TRUE))>0){
    bed1 <- bed[bed[,7]!=bed[,8],]  
    u1 <- GRanges(bed1[,1], IRanges(bed1[,2]+1, bed1[,7]), strand=bed1[,6], transcript_id=bed1[,4])
    u1 <- u1[start(u1)!=(end(u1)+1)]
    u2 <- GRanges(bed1[,1], IRanges(bed1[,8]+1, bed1[,3]), strand=bed1[,6], transcript_id=bed1[,4])
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
