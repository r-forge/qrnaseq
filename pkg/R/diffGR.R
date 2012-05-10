diffGR <- function(gr1, gr2, by="transcript_id"){
  idx <- match(by, colnames(gr1@elementMetadata))
  GR1 <- GRanges(paste(seqnames(gr1), gr1@elementMetadata[,idx], sep="|"), ranges(gr1), strand(gr1))
  GR2 <- GRanges(paste(seqnames(gr2), gr2@elementMetadata[,idx], sep="|"), ranges(gr2), strand(gr2))
  GRd <- suppressWarnings(setdiff(GR1, GR2))
  chrid <- do.call("rbind", strsplit(as.character(seqnames(GRd)), split="\\|"))
  GRD <- GRanges(seqnames=chrid[,1], ranges=ranges(GRd), strand(GRd), chrid[,2])
  colnames(GRD@elementMetadata) <- by
  return(GRD)
}
