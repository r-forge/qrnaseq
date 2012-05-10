\name{read.GenePred}
\alias{read.GenePred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  read gene annotation in GenePred format
}
\description{
  The function read GenePred format annotation data into GRange objects by
  different features. It can split the regions into different features,
  such as exon, CDS, intron and utr.
}
\usage{
read.GenePred(file, geneid = 12, feature = c("exon", "CDS", "intron",
"utr"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{The gene annotation file in GenePred format. The 9
    additional fields are required.}
  \item{geneid}{If not NULL, it is the column number of the input
    file.}
  \item{feature}{A character vector that are used to classify the gene
    regions. It can be the combination of "exon", "CDS", "intron" and
    "utr".}
  \item{...}{More options for read.table to load the GenePred file.}
}
\details{
  More details about the GenePred format can be found at:
  http://genome.ucsc.edu/FAQ/FAQformat#format9
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
  \code{\link{read.GTF}}, \code{\link{read.bed}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bed }% __ONLY ONE__ keyword per line
