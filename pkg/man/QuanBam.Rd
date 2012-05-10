\name{QuanBam}
\alias{QuanBam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Quantify aligned Bam file
}
\description{
  The function quantify a aligned Bam file from RNA-Seq with
  user-specified genome regions. The reads from bam file will be
  filtered by the regions. The output can be counts or maxdepths for
  each region.
}
\usage{
QuanBam(bam, GR, method = "both", Reduce="union", ovtype = "within", revstrand = FALSE, ingene = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam}{The aligned bam file from RNA-Seq experiment.}
  \item{GR}{\code{GRange} object with gene annotation.}
  \item{method}{The method to quantify reads. It can be count, maxdepth
    or both.}
  \item{Reduce}{The methods to reduce the overlapped regions. The
    regions from different transcripts but the same gene can be
    overlapped. If it is "union", the overlapped regions from the same
    gene will be united together. If it is "intersection", the
    returned object will be the most shared regions among the
    transcripts in the same gene. If it is NULL, all the regions will be
    quantified separately.}
  \item{ovtype}{The overlap type that reads map to genes. Only the reads
    mapping with the type can be used in the quantification. By default,
    the reads that aligned within the input regions is used. The other
    reads are filtered out. See more details in "type" section of
    help(findOverlaps).}
  \item{revstrand}{For paired reads. If TRUE, the paired reads that are
    in the same direction will be filtered out. If FALSE, the direction
    of reads will be ignored.}
  \item{ingene}{For paired reads. If TRUE, the paired reads that are not
    in the same gene will be filtered out. If FALSE, all the paired
    reads will be used.}
  \item{\dots}{More options when read bam file. See more details in
    help(readBamGappedAlignments)}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  The results return a GRanges object. If reduce option is used, the
  return GRanges object is reduced ranges by gene_id. Columns "counts"
  or/and "maxdepths" will be added to the elementMetadata of the GRanges
  object.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Qiang Hu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Quantification }
\keyword{ bam }% __ONLY ONE__ keyword per line
