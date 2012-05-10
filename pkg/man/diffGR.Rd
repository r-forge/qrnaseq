\name{diffGR}
\alias{diffGR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Perform the difference of two GRanges objects with gene annotation
}
\description{
  The function discard the overlapped regions of two GRanges objects
  from the first objects for the same gene.
}
\usage{
diffGR(gr1, gr2, by = "transcript_id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr1, gr2}{\code{GRanges} objects with annotation in their
    elementMetatdata.}
  \item{by}{The gene annotation where to perform the difference.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  The result returns a \code{GRanges} objects which is within gr1 but not
  overlapped with gr2 in each annotated level.
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
  \code{\link{setdiff}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Quantification }
\keyword{ setdiff }% __ONLY ONE__ keyword per line
