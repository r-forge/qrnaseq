\name{QuanGR}
\alias{QuanGR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summarize abundances to gene level
}
\description{
  The function summarizes the quantification results from the function QuanBam to gene
  level.
}
\usage{
QuanGR(GR, method = "sum", Quan = "counts", by = "gene_id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GR}{\code{GRange} object with quantification results and gene
    annotation.}
  \item{method}{The method to summarize abundances. It can be RPKM, sum,
    mean, median, max and min.}
  \item{Quan}{The quantification results from the function
    \code{QuanBam}. It can be counts or mdepths.}
  \item{by}{The gene annotation. The quantification results are
    summarized by this annotation level.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  The results return a abundance matrix in annotation level. For
  example, if the input abundances are summarized by gene_id level, the
  results will be a matrix with quantification result for each gene.
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
  \code{\link{QuanBam}}
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
\keyword{ RPKM }% __ONLY ONE__ keyword per line
