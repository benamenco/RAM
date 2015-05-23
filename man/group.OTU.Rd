\name{group.OTU}
\alias{group.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Distribution of OTUs
}
\description{
This function plot the distributon of otus in each level of a 
given metadata variable. The plot can be boxplot of barplot. 
The boxplot shows the range of relative abundance of a given 
otuID in each level of metadata category.  The barplot shows the 
relative abundance of the total counts of a given otuID in 
each level of metadata category. 
}
\usage{
group.OTU(otu, rank="g", otuIDs="", meta, meta.factor="", 
         boxplot=TRUE, main="", file=NULL, ext=NULL, 
         height=8, width=16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu}{
the OTU table to be analyzed.
}
  \item{rank}{
optional. if \code{NULL}, the lca of the OTUs will be retrieved; 
otherwise, the user should provide a valid taxonomic rank that an
otu being classified to (see ?RAM.rank.formatting for formatting 
details).
}
  \item{otuIDs}{
an vector of otuIDs in the OTU table 
}
  \item{meta}{
the metadata table to be analyzed.
}
  \item{meta.factor}{
the metadata qualitative variable 
}
  \item{boxplot}{
logical. If \code{TRUE}, generate boxplot; 
otherwise generate barplot.
}
  \item{main}{
title of the plot.
}
  \item{file}{
the file path where the image should be created (see ?RAM.plotting).
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
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{group.OTU} returns boxplot or barplot for the distribution 
of a list of otuIDs.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
\note{
The OTUs are determined to be absent/present using the "pa" method 
from the function \code{\link{decostand}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ggplot}}
}
\examples{
data(ITS1, meta)

# otuIDs
otuIDs=rownames(ITS1)[1:10]

# names(meta)
theme <- RAM.color()
group.OTU(otu=ITS1, rank="g", otuIDs=otuIDs, 
        meta=meta, meta.factor="City", boxplot=TRUE,  
        file=NULL, ext=NULL) + theme

\dontrun{
group.OTU(otu=ITS1, rank="g", otuIDs=otuIDs, 
               meta=meta, meta.factor="City", boxplot=FALSE, 
               file=NULL) + theme
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
