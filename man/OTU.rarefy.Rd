\name{OTU.rarefy}
\alias{OTU.rarefy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Rarefied OTU Tables
}
\description{
This function output rarefied OTU tables using \code{\link{rrarefy}}.
This function may take long time for large dataset, e.g. over 100k
otus x 45 samples.
}
\usage{
OTU.rarefy(data, sample=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of otu tables. See also \code{\link{RAM.input.formatting}}.
}
  \item{sample}{
an integer represent the sampling size.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function returns a list of rarefied otu tables.
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
\note{
See also \code{\link{rrarefy}} 
}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link{RAM.input.formatting}}.
}
\examples{
\dontrun{
data(ITS1, ITS2)
otus.rf <- OTU.rarefy(data=list(ITS1=ITS1, ITS2=ITS2), 
                      sample=NULL)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
