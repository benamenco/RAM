\name{shared.Taxa}
\alias{shared.Taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of Shared Taxa Across ALL Subjects
}
\description{
This function consumes OTU tables or a taxonomy matrices and 
returns a list summarizing information about the presence of 
the taxa in that table at a given taxonomic rank.
}
\usage{
shared.Taxa(data, is.OTU=TRUE, rank="g")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables or taxonomy abundance matrices to be analyzed.
}
  \item{is.OTU	}{
whether or not the input data are otu tables
}
  \item{rank}{
the taxonomic rank to be investigated
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{shared.Taxa} returns a list containing the information 
calculated. The names associated with the list describe what that 
number represents; i.e. "#_of_families_in_all_samples" shows how 
many taxa at the family level were found to be present in all samples. 
The last item in the list is a character vector, containing the 
taxon names of which were present in all samples.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
\note{
The taxa are determined to be absent/present using the "pa" method 
from the function \code{\link{decostand}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{decostand}}
}
\examples{
data(ITS1)
shared.Taxa(data=list(ITS1=ITS1))
\dontrun{
g1 <- tax.abund(ITS1, rank="g", drop.unclassified=TRUE)
shared.Taxa(data=list(genus_ITS1=g1), rank="g", is.OTU=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
