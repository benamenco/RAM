\name{OTU.diversity}
\alias{OTU.diversity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize Diversity Indices for OTU Tables
}
\description{
These functions calculate diversity indices for all samples and 
append outputs as new columns to metadata table. 
}
\usage{
OTU.diversity(data, meta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables. 
}
  \item{meta}{
the metadata to append the outputs.
}
}
\details{
This function summarize the following diversity indices: specnumber, 
shannon, simpson, invsimpson, true diversity, evenness, chao and 
ACE indices, for a given otu table. 
See \code{\link{diversity.indices}})
}
\value{
This function return vectors of diversity indices for each sample, 
which are appended to a given metadata table. 
}
\note{
Credit goes to package \code{vegan} for the partial argument matching 
(see References), and for the calculation of all diversity indices 
except for true diversity and evenness.
}
\references{
Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre, 
Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson,
  Peter Solymos, M. Henry H. Stevens and Helene Wagner (2013). 
  vegan: Community Ecology Package. R package version
  2.0-10. http://CRAN.R-project.org/package=vegan
  
Diversity index. (2014, May 7). In Wikipedia, The Free Encyclopedia. 
Retrieved 14:57, May 28, 2014, from
http://en.wikipedia.org/w/index.php?title=Diversity_index&oldid=607510424
  
Blackwood, C. B., Hudleston, D., Zak, D. R., & Buyer, J. S. (2007). 
Interpreting ecological OTU.diversity applied to terminal restriction 
fragment length polymorphism data: insights from simulated microbial 
communities. Applied and Environmental Microbiology, 73(16), 5276-5283.
}
\author{
Wen Chen.
}
%\note{
%}
\examples{
data(ITS1, ITS2, meta)
data=list(ITS1=ITS1, ITS2=ITS2)
\dontrun{
meta.diversity<-OTU.diversity(data=data, meta=meta) 
}                              
}
\keyword{manip}
