
QuanGR <- function(GR, method="sum", Quan="counts", by="gene_id"){
  METHODS <- c("RPKM", "sum", "mean", "median", "max", "min")
  idx <- pmatch(method, METHODS)
  qidx <- pmatch(Quan, colnames(GR@elementMetadata))
  gidx <- pmatch(by, colnames(GR@elementMetadata))
  if(is.na(qidx)){
    stop(paste(Quan, "is not available"))
  }
  RPKM <- function(GR){
    N <- sum(GR@elementMetadata[,qidx])
    C <- tapply(GR@elementMetadata[,qidx], GR@elementMetadata[,gidx], sum)
    L <- tapply(width(GR), GR@elementMetadata[,gidx], sum)
    rpkm <- (10^9*C)/(N*L)
    return(rpkm)
  }

  if(is.na(idx)){
    stop("invalid method")
  }else if(idx==1){
    quan <- RPKM(GR)
  }else{
    quan <- tapply(GR@elementMetadata[,qidx], GR@elementMetadata[,gidx], METHODS[idx]) 
  }
  quan <- cbind(quan)
  colnames(quan) <- paste(METHODS[idx], "_", Quan, sep="")
  return(cbind(quan))
}
