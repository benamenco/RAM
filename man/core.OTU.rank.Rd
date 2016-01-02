\name{core.OTU.rank}
\alias{core.OTU.rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary Of Core OTUs
}
\description{
This function returns a list showing otus that present in a
pre-defined percent of samples in each level of a given metadata
category.
}
\usage{
core.OTU.rank(data, rank="g", drop.unclassified=TRUE,
                meta, meta.factor="", percent=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables to be analyzed.
See also \code{\link{RAM.input.formatting}}.
}
  \item{rank}{
the taxonomic rank(s) of otu classification
(see ?RAM.rank.formatting for
formatting details).
}
  \item{drop.unclassified}{
logical, whether or not exclude unclassified groups.
}
  \item{meta}{
the metadata table to be analyzed.
}
  \item{meta.factor}{
the metadata qualitative variable
}
  \item{percent}{
the percent of samples in each level of the given metadata
variable
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{core.OTU.rank} returns a list containing otus that present
in a pre-defined percent of samples in each level of a given 
metadata category. The outputs describe the following 
information for each level of a given metadata variable: 
1) core otuID; 2) taxa the core otus assigned to; and 
3) percent of core otus sequences vs. total sequences in each 
group. The last item in the list show the same information of 
otus that in all levels.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
\note{
The taxon groups are determined to be absent/present using the 
"pa" method from the function \code{\link{decostand}}.
}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~
\seealso{
\code{\link{decostand}}
}
\examples{
data(ITS1, meta)
\dontrun{
core <- core.OTU.rank(data=list(ITS1=ITS1), rank="g", meta=meta,
         meta.factor="City", percent=0.90)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
