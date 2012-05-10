\name{read.GTF}
\alias{read.GTF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  read gene annotation in GTF format
}
\description{
  The function read GTF format annotation data into GRange objects by
  different features. The input GTF data can be downloaded from UCSC,
  Ensembl and so on. The function can split the regions into
  different features, such as exon, CDS, intron and utr.
}
\usage{
read.GTF(file, feature = c("exon", "CDS", "intron", "utr"),
attributes = c("gene_id", "transcript_id", "gene_name",
"transcript_name", "protein_id")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{The gene annotation file in GTF format.}
  \item{feature}{A character vector that are used to classify the gene
    regions. It can be the combination of "exon", "CDS", "intron" and
    "utr".}
  \item{attributes}{The element names of the attribute column in the GTF
  data.}
}
\details{
  More details about the GTF format can be found at:
  http://mblab.wustl.edu/GTF22.html
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  The results return a list of GRanges object according to the specified
  features. By default, the features "exon", "CDS", "intron", "utr" are
  specified. The results will contain a list of regions that contains
  exon, coding sites, intron, 3 end UTR and 5 end UTR from the input GTF
  file.
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
  \code{\link{read.bed}}, \code{\link{read.GenePred}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GTF }% __ONLY ONE__ keyword per line
