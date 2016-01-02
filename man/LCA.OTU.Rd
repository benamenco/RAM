\name{LCA.OTU}
\alias{LCA.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lowest Common Ancestor (LCA) OF EACH OTU
}
\description{
This function consumes an OTU table and extract the LCA (lowest
common ancestor) that each otu assigned to.
See also \code{\link{tax.split}}.
}
\usage{
LCA.OTU(otu, strip.format=FALSE, drop=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu}{
the OTU table to be used.
}
  \item{strip.format}{
logical. Whether or not to remove the prefix of the taxonomy
assignment at each rank.
see
}
  \item{drop}{
logical. Whether or not drop taxonomic columns other than LCA.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function return a data frame same as the input OTU table,
except the last column is the LCA of each otu, not the lineage.
The taxonomy column can be kept, by using \code{drop}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
\note{
\code{tax.split} returns the same OTU table with classification
at a given taxonomic rank, which can be missing if an otu was 
not classified a that that taxonomic level. 
\code{\link{LCA.OTU}}, guaranteed that all OTUs will be 
preserved in the returned data table and provide the LCA for 
each OTU, although only higher taxonomic ranks were available.
}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~
\seealso{
\code{\link{tax.split}}
}
\examples{
data(ITS1)
\dontrun{
# compare the following 2 commands:
# keep the rank prefix of the LCA column
ITS1.LCA <- LCA.OTU(ITS1, strip.format=TRUE, drop=TRUE)
# remove the rank prefix of the LCA column
ITS1.LCA <- LCA.OTU(ITS1, strip.format=FALSE, drop=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
