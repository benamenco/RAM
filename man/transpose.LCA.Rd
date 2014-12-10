\name{transpose.LCA}
\alias{transpose.LCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transpose OTU Tables With LCA Annotation For Each OTU
}
\description{
Similar to \code{\link[RAM]{transpose.OTU}}, but each OTU is 
annotated by the lowest common ancestor it was assigned to.
}
\usage{
transpose.LCA(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The OTU tables to be transposed.
See also \code{\link{RAM.input.formatting}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a transposed OTU table, but the colname is formatted as:
LCA_otuID.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, ITS2)
\dontrun{
lca.t <- transpose.LCA(data=list(ITS1=ITS1, ITS2=ITS2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
