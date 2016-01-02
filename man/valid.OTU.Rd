\name{valid.OTU}
\alias{valid.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validate an OTU Table
}
\description{
This function consumes one or two OTU tables and checks if they
are formatted properly and contain valid data.
}
\usage{
valid.OTU(otu1, otu2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to check.
}
  \item{otu2}{
the second OTU table to check.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
If the table is not valid, an error will be raised with a
description explaining the problem. If the table is valid, NULL
will be returned invisibly.
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
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, ITS2)
valid.OTU(ITS1)
valid.OTU(ITS1, ITS2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{error}
