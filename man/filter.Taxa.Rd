\name{filter.Taxa}
\alias{filter.Taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Taxonomic Abundance Matrix by Total Counts Or
Maximum Relative Abundance
}
\description{
This function filter taxa group by counts or relative abundance.
If filter by counts, taxa having total counts more than a 
threshhold will be kept. If filter by relative abundance, taxa
with the maximum relative abundance greater than a threshhold 
in at least one subject will be kept.
}
\usage{
filter.Taxa(taxa, drop.unclassified=TRUE,
            percent=NULL, number=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxa}{
the taxonomy abundance matrix: sample x species data frame.
See also \code{\link{tax.abund}}
}
  \item{drop.unclassified}{
logical, whether or not remove unclassified groups.
See also \code{\link{tax.abund}}
}
  \item{percent}{
a floating point greater than 0 and less or equals to 1.
Filter Taxa table by relative abundance.
}
  \item{number}{
an integer. FilterTaxa table by total sequence counts.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a data frame with taxa 
met the filter requirement only.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
% \note{
%% ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
%\seealso{
% ~~further links~~
%}
\examples{
data(ITS1)
g1 <- tax.abund(ITS1, rank="g", drop.unclassified=TRUE)
taxa.fil <- filter.Taxa(g1, percent=0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
