\name{transpose.OTU}
\alias{transpose.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Take the Transpose of an OTU Table
}
\description{
Returns the transpose of the given OTU table, excluding the last 
column (which should contain taxonomic information).
}
\usage{
transpose.OTU(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The OTU table to be transposed.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a data frame with rows equal to the columns of the original 
OTU, and columns equal to the rows of the original OTU. 
(Excluding the taxonomy column).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1)

ITS1.t <- transpose.OTU(ITS1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
