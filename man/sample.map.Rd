\name{sample.map}
\alias{sample.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot The Geographic Location of Samples
}
\description{
This function plot the number of samples collected from different 
locations that are DISTANT from each other, e.g. samples that 
collected from distant cities. This function is similar
but not the same as \code{\link{sample.locations}} and 
\code{\link{sample.sites}}.
The plot will also show the sample size of each location.
}
\usage{
sample.map(meta, siteID="City", maptype="roadmap",  
           lat="Latitude", lon="Longitude", zoom=3, 
           file=NULL, ext=NULL, width=10, height=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta}{
the OTU table to be used.
}
  \item{siteID}{
IDs of sampling sites for each unique pair of longitude and 
latitude. 
}
  \item{maptype}{
maptype to use, see also \code{\link{get_map}}.
}
  \item{lat}{
latitude of each sampling location
}
  \item{lon}{
longitude of each sampling location
}
  \item{zoom}{
map zoom, an integer from 3 (continent) to 21 (building). 
see also \code{\link{get_map}}.
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
Please note that this function is getting map information using 
either the Google Maps API or the OpenStreetMap API, and your 
usage is subject to the terms of those APIs. 
}

\note{
If you are getting a 403/503 error, that likely means that the 
current map provider is currently unavailable. Try calling the 
function again in a few hours.
If you get a warning message of the form "\code{Removed X rows 
containing missing values (geom_point).}", this means that the 
current zoom level is too high to display some or all of the points. 
Try using a lower value for\code{zoom}.
}
%\value{
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sample.locations}}, \code{\link{sample.sites}}
}
\examples{
data(meta)

\dontrun{
sample.map(meta=meta, zoom=8)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
