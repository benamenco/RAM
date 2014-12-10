\name{percent.classified}
\alias{percent.classified}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Percent of OTUs Classified at a Given Taxonomic Rank
}
\description{
This function consumes an OTU table, and a vector containing 
taxnomic ranks, then returns what percent of OTUs in the given 
table are classified at each taxonomic rank.
}
\usage{
percent.classified(data, ranks=c("f","g"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables to be processed. 
See also \code{\link{RAM.input.formatting}}
}
  \item{ranks}{
a vector containing the taxonomic ranks you are interested in 
(see ?RAM.rank.formatting for formatting details).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of numeric vectors, containing the result for each taxonomic rank.
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
data(ITS1, ITS2)
data <- list(ITS1=ITS1, ITS2=ITS2)
# find what percent of OTUs classified at family and genus levels
percent.classified(data=data, ranks=c("f","g"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
