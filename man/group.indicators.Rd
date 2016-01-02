\name{group.indicators}
\alias{group.indicators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Indicator Taxon Groups for Metadata Trends
}
\description{
This function conumes one or two OTU tables, along with a 
metadata facotr, and creates a barplot showing the relative 
abundance of all groups which are statistical indicators of 
that factor.  
}
\usage{
group.indicators(data, is.OTU=TRUE, meta, factor, rank,
                thresholds = c(A = 0.85,
                               B = 0.8,
                               stat = 0.8,
                               p.value = 0.05),
                 cex.x=NULL, file = NULL, ext = NULL, 
                 height = 12, width = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTUs or taxonomic abundance matrices.
see also \code{\link{RAM.input.formatting}}
}
  \item{is.OTU}{
logical. Whether the input data are OTU tables.
}
  \item{meta}{
the metadata table to be used.
}
  \item{factor}{
a named character vector (of length 1) giving the name of the
column in \code{meta} to be used when performing the analysis
(see \code{\link{RAM.factors}}).
}
  \item{rank}{
the taxonomic rank to use
(see ?RAM.rank.formatting for formatting details).
if \code{rank} is \code{NULL}, will use otus as indicators 
which are annotated by the lca assigned to the otus, otherwise 
will use taxon names as indicators at the given taxonomic rank.
}
  \item{thresholds}{
a character vector of length 4 specifying the thresholds for 
the A, B, stat, and p values (see Details).
}
  \item{cex.x}{
optional. The size of x axis names.
}
  \item{file}{
the file path where the image should be created
(see ?RAM.plotting).
}
  \item{ext}{
the file type to be used; one of \code{"pdf"}, \code{"png"},
\code{"tiff"}, \code{"bmp"}, \code{"jpg"}, or \code{"svg"}.
}
  \item{height}{
the height of the image to be created (in inches).
}
  \item{width}{
the width of the image to be created (in inches).
}
}
\details{
This function uses \code{indicspecies::multipatt} to determine
the indicators. After this analysis is performed, there will
likely be some species determined to be 'significant,' but to
varying degrees. To control how many groups are selected, you
can adjust the \code{thresholds} argument. It consists of four
components:
(quotations taken from \code{vignette("indicspeciesTutorial")},
see References):
\describe{
  \item{A}{the \emph{specificity} of the indicator; 'the
  probability that the surveyed site belongs to the target 
  site   group given the fact that the species has been found'}
  \item{B}{the \emph{fidelity} of the indicator; 'the 
  probability of finding the species in sites belonging to the 
  site group'}
  \item{stat}{the association strength for the combinations.}
  \item{p.value}{'the degree of statistical significance of
  the association'}
}
Only the species with \code{A}, \code{B}, and \code{stat} 
values above, and \code{p.value} below those given in 
\code{thresholds} will be kept.
}
\value{
This function returns a stacked barplot and a vector of 
identified indicators, including the ones in the plot and the 
ones being excluded from the plot.
}
\references{
De Caceres, M., Legendre, P. (2009). Associations between 
species and groups of sites: indices and statistical inference.
Ecology, URL http://sites.google.com/site/miqueldecaceres/
}
\author{
Wen Chen and Joshua Simpson.
}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~
\seealso{
\code{\link[indicspecies]{multipatt}}
}
\examples{
data(ITS1, ITS2, meta)
\dontrun{
# inputs are OTU tables
group.indicators(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE, 
                 meta, factor = c(Province="Province"),
                 rank="g")
group.indicators(data=list(ITS1=ITS1), is.OTU=TRUE, meta,
                     factor = c(Province="Province"),
                     rank=NULL)
group.indicators(data=list(ITS1=ITS1), is.OTU=TRUE, meta,
                     factor = c(Province="Province"),
                     rank=NULL)
# inputs are taxonomic abundance matrices
g1 <- tax.abund(ITS1, rank="g")
g2 <- tax.abund(ITS2, rank="g")
group.indicators(data=list(g1=g1, g2=g2), is.OTU=FALSE, meta,
                     factor = c(Province="Province"),
                     rank="g")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
