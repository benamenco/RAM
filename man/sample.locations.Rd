\name{sample.locations}
\alias{sample.locations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Geographic Location of Samples
}
\description{
This function consumes an OTU table, along with its associated metadata, and 
plots all the samples from that data as points on a map. The size of a point 
indicates the number of counts collected from that sample, while the colour of 
the point (optionally) shows a metadata factor for that sample.
}
\usage{
sample.locations(otu1, otu2=NULL, meta, factor = NULL, zoom = 5, 
                 source = "google", labels = c("ITS1", "ITS2"),
                 lat.col = "Latitude", long.col = "Longitude",
                 file = NULL, ext = NULL, height = 10, width = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the OTU table to be used.
}
  \item{otu2}{
the (optional) second OTU table to be used.  
}
  \item{meta}{
the metadata table to be used.
}
  \item{factor}{
(optional) a named character vector of length one specifying a column from the
metadata table to be used to colour the points.
}
  \item{zoom}{
a positive integer in the range 3-21 (if \code{source == "google"}) or 3-18 (if
\code{source == "osm"}) specifying the zoom for the map (low number means zoomed
out).
}
  \item{source}{
the source to be used for the map; either "\code{google}" or "\code{osm}". 
}
\item{labels}{
a character vector giving one label per OTU.
}
  \item{lat.col}{
the name of the column in \code{meta} containing the latitude information.
}
  \item{long.col}{
the name of the column in \code{meta} containing the longitude information.
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
Please note that this function is getting map information using either the 
Google Maps API or the OpenStreetMap API, and your usage is subject to the
terms of those APIs. 
}

\note{
If you are getting a 403/503 error, that likely means that the current map
provider is currently unavailable. This can be for a variety of reasons:
if \code{source == "google"}, you have likely maxed out your API call limit
(this can be due to multiple users sharing an IP address; contact your 
system administrator for further details). If \code{source == "osm"}, the 
server is likely under heavy load and unable to process your request. You
can check the server load 
\href{http://munin.openstreetmap.org/openstreetmap/yevaud.openstreetmap/load.html}{online}. 
In either case, the issue will likely resolve itself. Try
calling the function again in a few hours.

If you get a warning message of the form "\code{Removed X rows containing 
missing values (geom_point).}", this means that the current zoom level is too 
high to display some or all of the points. Try using a lower value for
\code{zoom}.
}
%\value{
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{RAM.factors}}
}
\examples{
data(ITS1, meta)

\dontrun{
sample.locations(ITS1, meta, factor=c(Crop="Crop"))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
