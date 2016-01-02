\name{group.rich}
\alias{group.rich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barplot Of Richness For Each Level Of A Given Metadata Variable
}
\description{
This function first use \code{\link{specpool}} to estimate the
extrapolated species richness in a species pool (levels of 
metadata variable), and the number of unobserved species, then 
do a barplot.
}
\usage{
group.rich(otu, meta, factor, file=NULL, ext=NULL, 
           width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu}{
an OTU table.
}
  \item{meta}{
the metadata table to be used (must have same samples as
\code{data}.
}
  \item{factor}{
a character string. Must be one of the metadata variables.
}
  \item{file}{
optional. Filename that the plot to be saved to.
}
  \item{ext}{
optional. Filename extension, type of image to be saved.
}
  \item{width}{
an integer. Filter OTU table by counts.
}
  \item{height}{
an integer. Filter OTU table by counts.
}
}
%\details{
%%  Maybe explain full/restricted models
%}
\value{
This function returns a barplot of species richness for a given
metadata variable.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
%\note{
%%  ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link[vegan]{specpool}}, \code{\link[vegan]{specpool}}
}
\examples{
data(ITS1, meta)
group.rich(ITS1, meta, "Crop")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
