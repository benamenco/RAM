\name{group.indicators}
\alias{group.indicators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Indicator Taxon Groups for Metadata Trends
}
\description{
This function conumes one or two OTU tables, along with a metadata facotr, and 
creates a barplot showing the relative abundance of all groups which are 
statistical indicators of that factor.
}
\usage{
group.indicators(otu1, otu2 = NULL, meta, factor, rank, 
                 thresholds = c(A = 0.85, B = 0.8, stat = 0.8, p.value = 0.05), 
                 labels = c("ITS1", "ITS2"), 
                 file = NULL, ext = NULL, height = 12, width = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be used.
}
  \item{otu2}{
the second OTU table to be used.
}
  \item{meta}{
the metadata table to be used.
}
  \item{factor}{
a named character vector (of length 1) giving the name of the column in
\code{meta} to be used when performing the analysis (see
\code{\link{RAM.factors}}).
}
  \item{rank}{
the taxonomic rank to use (see ?RAM.rank.formatting for formatting details).
}
  \item{thresholds}{
a character vector of length 4 specifying the thresholds for the A, B, stat, and
p values (see Details).
}
  \item{labels}{
a character vector giving the labels for the OTU tables.
}
  \item{file}{
the file path where the image should be created (see ?RAM.plotting).
}
  \item{ext}{
the file type to be used; one of \code{"pdf"}, \code{"png"}, \code{"tiff"}, 
\code{"bmp"}, \code{"jpg"}, or \code{"svg"}.
}
  \item{height}{
the height of the image to be created (in inches).
}
  \item{width}{
the width of the image to be created (in inches).
}
}
\details{
This function uses \code{indicspecies::multipatt} to determine the indicators.
After this analysis is performed, there will likely be some species determined
to be 'significant,' but to varying degrees. To control how many groups are
selected, you can adjust the \code{thresholds} argument. It consists of four
components (quotations taken from \code{vignette("indicspeciesTutorial")}, see
References):

\describe{
  \item{A}{the \emph{specificity} of the indicator; "the probability that the
surveyed site belongs to the target site group given the fact that the species
has been found."}
  \item{B}{the \emph{fidelity} of the indicator; "the probability of finding the
species in sites belonging to the site group."}
  \item{stat}{the association strength for the combinations.}
  \item{p.value}{"the degree of statistical significance of the association."}
}

Only the species with \code{A}, \code{B}, and \code{stat} values above, and 
\code{p.value} below those given in \code{thresholds} will be kept.
}
%\value{

%}
\references{
De Caceres, M., Legendre, P. (2009). Associations between species and groups of
sites: indices and statistical inference. Ecology, URL
http://sites.google.com/site/miqueldecaceres/

}
\author{
Wen Chen and Joshua Simpson.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multipatt}}
}
\examples{
data(ITS1, ITS2, meta)

group.indicators(ITS1, ITS2, meta, factor = c(Province="Province"), rank="g")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
