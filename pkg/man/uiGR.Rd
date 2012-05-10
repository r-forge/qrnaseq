\name{uiGR}
\alias{uiGR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reduce GRange object in gene level
}
\description{
  The function reduce overlapped GRanges objects with union or
  intersection method in gene level. Usually, distinct transcripts of a
  gene share several the same exon regions. The function can reduce
  these overlapped exon into union or intersection regions.
}
\usage{
uiGR(GR, Reduce = "union")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GR}{\code{GRange} object with gene_id annotation.}
  \item{Reduce}{The methods to reduce the overlapped regions. The
    regions from different transcripts but the same gene can be
    overlapped. If it is "union", the overlapped regions from the same
    gene will be united together. If it is "intersection", the
    returned object will be the most shared regions among the
    transcripts in the same gene.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  The results return a GRanges object with gene_id annotation. If
  "union" is used, the results are the united regoins within each
  gene. If the "intersection" is used, the regions are the most shared
  regions in each gene.
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
\keyword{ union }
\keyword{ intersection }% __ONLY ONE__ keyword per line
