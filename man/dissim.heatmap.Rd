\name{dissim.heatmap}
\alias{dissim.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Distance Matrix Heatmap for OTU Samples
}
\description{
This function consumes an OTU table, along with some optional 
parameters, and creates a heatmap showing the distance matrix 
for the samples of the given OTU table.
}
\usage{
dissim.heatmap(data, meta=NULL, row.factor=NULL, col.factor=NULL, 
               stand.method="chi.square", dissim.method="euclidean",
               file=NULL, ext=NULL, height=8, width=9, 
               leg.x=-0.05, leg.y=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be used.
}
  \item{meta}{
the metadata table to be used.  
}
  \item{row.factor}{
a factor from the metadata to show along the rows of the heatmap 
(see Details below).
}
  \item{col.factor}{
a factor from the metadata to show along the columns of the heatmap 
(see Details below).
}
  \item{stand.method}{
a method used to standardize the OTU table. One of "\code{total}", 
"\code{max}", "\code{freq}", "\code{normalize}", "\code{range}", 
"\code{standardize}", "\code{pa}", "\code{chi.square}", 
"\code{hellinger}" or "\code{log}" (see ?decostand).
}
  \item{dissim.method}{
the dissimilarity index to be used; one of "\code{manhattan}",
"\code{euclidean}", "\code{canberra}", "\code{bray}", "\code{kulczynski}",
"\code{jaccard}", "\code{gower}", "\code{altGower}","\code{morisita}",
"\code{horn}", "\code{mountford}", "\code{raup}",
"\code{binomial}","\code{chao}", or "\code{cao}" (see ?vegdist).
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
  \item{leg.x}{
how far the legend should be inset, on the x axis.
}
  \item{leg.y}{
how far the legend should be inset, on the y axis.
}
}
\details{
Both \code{row.factor} and \code{col.factor} should be named 
character vectors specifying the names of the rows to be used from 
meta (see \code{\link{RAM.factors}}). They should also be factors; 
if they are not, a warning is raised and they are coerced to factors 
(see \code{\link{factor}}). A warning is also raised when a factor 
has more than 8 levels, as that is the most colours the current 
palettes support. The factor must also correspond to the OTU table; 
i.e. they should have the same samples. If not, an error is raised.
}

\note{
This function creates the heatmap using the \code{heatmap.2} 
function from the \code{gplots} package. That function calls 
\code{layout} to set up the plotting environment, which currently 
prevents plotting two data sets side by side, or to automatically 
place the (metadata) legend. Unfortunately, this means that the
\code{leg.x} and \code{leg.y} parameters must be used to adjust 
the legend by trial and error. It is possible to move the legend 
outside of the plotting area; if not legend appears, try with 
small \code{leg.x} and \code{leg.y} values.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{decostand}}, \code{\link{vegdist}}, \code{\link{factor}},
\code{\link{heatmap.2}}, \code{\link{RAM.factors}}
}
\examples{
data(ITS1, meta)

# plot to the screen with one meta factor and standard 
# calculation methods 
dissim.heatmap(ITS1, meta, row.factor=c(Plot="Plots"))

\dontrun{
# plot the heatmap to a .tiff file using Hellinger standardization 
# and Manhattan distances
dissim.heatmap(ITS1, dissim.method=manhattan, 
               stand.method=hellinger,
               file="my/sample/path", ext="tiff")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
