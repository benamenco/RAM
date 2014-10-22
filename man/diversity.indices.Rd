\name{diversity.indices}
\alias{diversity.indices}
\alias{true.diversity}
\alias{evenness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate True Diversity and Evenness 
}
\description{
These functions calculate true diversity and evenness for all samples.
}
\usage{
true.diversity(otu1, otu2 = NULL, index = "simpson")
evenness(otu1, otu2 = NULL, index = "simpson")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the otu1 data to be processed.
}
  \item{otu2}{
the otu2 data to be processed.
}
  \item{index}{
the index to use for calculations; partial match to "\code{simpson}" or
"\code{shannon}".  
}
}
\details{
The \code{otu2} argument is optional; if it is not specified, the
 calculations will be done for only the \code{otu1} argument.

For the following sections, \eqn{S}{S} represents the number of 
species, \eqn{\lambda}{\lambda} represents the Simpson index, and 
\eqn{H'}{H'} represents the Shannon index.

The formulas for the true diversity of the indices are as follows:

\itemize{
  \item Simpson: \eqn{D_2 = frac{1}{\lambda}}{D_2 = 1 / \lambda}
  \item Shannon: \eqn{D_1 = \exp{H'}}{D_1 = exp(H')}
}

The formulas for the evenness of the indices are as follows:

\itemize{
  \item Simpson: \eqn{\frac{\frac{1}{\lambda}}{S}}{(1 / \lambda) / S}
  \item Shannon: \eqn{\frac{H'}{\ln{S}}}{H' / ln(S)}
}
}
\value{
Both functions return a numeric matrix, where the rows are the given 
OTUs, and the columns are the samples.
}
\note{
Credit goes to package \code{vegan} for the partial argument matching 
(see References).
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
Interpreting ecological diversity indices applied to terminal restriction 
fragment length polymorphism data: insights from simulated microbial 
communities. Applied and Environmental Microbiology, 73(16), 5276-5283.
}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%}
\examples{
data(ITS1, ITS2)

# true diversity, using default index (Simpson)
true.diversity(ITS1)

# true diversity for ITS1 and ITS2, using Shannon
true.diversity(ITS1, ITS2, index="shannon")

# default evenness (Simpson) for ITS1/ITS2
evenness(ITS1, ITS2)

# Shannon evenness
evenness(ITS1, index="shannon")
}
\keyword{manip}
