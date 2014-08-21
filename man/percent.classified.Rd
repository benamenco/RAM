\name{percent.classified}
\alias{percent.classified}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Percent of OTUs Classified at a Given Taxonomic Rank
}
\description{
This function consumes an OTU table, and a rank, then returns what percent of
OTUs in the given table are classified at that taxonomic rank.
}
\usage{
percent.classified(data, rank)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be processed.
}
  \item{rank}{
the taxonomic rank you are interested in (see ?RAM.rank.formatting for
formatting details).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector of length one, containing the result.
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

# find what percent of OTUs classified at genus level
percent.classified(ITS1, "g")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
