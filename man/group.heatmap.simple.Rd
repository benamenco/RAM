\name{group.heatmap.simple}
\alias{group.heatmap.simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a Heatmap Showing OTU Abundance by Taxonomic Classification
}
\description{
This function consumes an OTU table and a rank, as well as some
optional parameters, and creates a heatmap showing the abundance
of the OTUs at the given taxonomic rank for each sample.
}
\usage{
group.heatmap.simple(data, is.OTU=TRUE, meta=NULL, rank, 
                     row.factor=NULL, top=NULL, count=FALSE, 
                     drop.unclassified=FALSE,
                     dendro="none", file=NULL, ext=NULL,
                     width=9, height=8, leg.x=-0.08, leg.y=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be used.
}
  \item{is.OTU}{
logical. Whether or not the data is an OTU table.
}
  \item{meta}{
the metadata table to be used.
}
  \item{rank}{
the taxonomic rank to use (see ?RAM.rank.formatting for 
formatting details).
}
  \item{row.factor}{
a factor from the metadata to show along the rows of the heatmap.
(see Details below).
}
  \item{dendro}{
a character vector specifying on which axes (if any) a 
dendrogram should be plotted. Must be one of "\code{none}", 
"\code{both}", "\code{column}", or "\code{row}".
}
  \item{top}{
the number of groups to select, starting with the most abundant.
If \code{NULL}, all are selected.
}
  \item{count}{
logical. Should the actual count of each OTU be shown, or should
the relative abundances be shown?
}
  \item{drop.unclassified}{
logical. Should OTUs labelled "unclassified" or missing 
classification at the given taxonomic rank be excluded?
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
  \item{leg.x}{
how far the legend should be inset, on the x axis.
}
  \item{leg.y}{
how far the legend should be inset, on the y axis.
}
}
\details{
\code{row.factor} should be a named character vector specifying 
the name of the row to be used from \code{meta} 
(see \code{\link{RAM.factors}}).

It should also be a factor; if it is not, a warning is raised 
and it is coerced to a factor (see \code{\link{factor}}). 
A warning is also raised when a factor has more than 8 levels, 
as that is the most colours the current palettes support. The 
factor must also correspond to the OTU table; i.e. they should 
have the same samples. If not, an error is raised.
}
\note{
This function creates the heatmap using the \code{heatmap.2} 
function from the \code{gplots} package. That function calls 
\code{layout} to set up the plotting environment, which 
currently prevents plotting two data sets side by side, or to 
automatically place the (metadata) legend. Unfortunately, this 
means that the \code{leg.x} and \code{leg.y} parameters must be 
used to adjust the legend by trial and error. It is possible to
move the legend outside of the plotting area; if no legend 
appears, try with small \code{leg.x} and \code{leg.y} values.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link{factor}}, \code{\link{heatmap.2}}, 
\code{\link{RAM.factors}}
}
\examples{
data(ITS1, meta)
# plot the abundance of all observed classes for each sample, 
# displaying it to the screen and adding a dendrogram on the 
# column and the Collector on the row
group.heatmap.simple(ITS1, is.OTU=TRUE, meta=meta,
                     row.factor=c(Crop="Crop"), dendro="row",
                     rank="g", top=10, drop.unclassified=TRUE,
                     leg.x=-0.06)
\dontrun{
# plot the genus for all OTUs, add a dendrogram to the row and 
# column, and save the plot in path.tiff
group.heatmap.simple(ITS1, is.OTU=TRUE, meta=meta, rank="genus",
                     dendro="both", file="my/file/path")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
