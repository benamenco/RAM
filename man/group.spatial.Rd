\name{group.spatial}
\alias{group.spatial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Spatial Collection Trends for Taxon Groups
}
\description{
This function consumes an OTU table and its associated metadata 
table, and uses that information to produce a choropleth (essentially 
a heatmap, but superimposed imposed on an actual, cartographic map) 
to show how many counts of each taxon group were detected in each 
Canadian province/territory.
}
\usage{
group.spatial(data, meta, date.col, province.col, rank, group, 
              breaks = "year", file = NULL, ext = NULL, 
              height = 8, width = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be used.
}
  \item{meta}{
the metadata table to be used.
}
  \item{date.col}{
a character vector specifying which column in metadata contains 
the date information (see \code{\link{RAM.dates}}).
}
  \item{province.col}{
a character vector specifying which column in metadata contains 
the province information (see Details).
}
  \item{rank}{
a character vector specifying the rank of the desired taxon groups. 
Note that all groups should come the same rank. 
(see \code{\link{RAM.rank.formatting}}).
}
  \item{group}{
a character vector giving the names of the groups to be plotted.
}
  \item{breaks}{
how many time segments should be plotted; see Details.
}
  \item{file}{
the file path where the image should be created. 
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
This function currently only supports Canadian data. The entries 
in meta$province.col should be specified as provinces, using the 
standard \href{http://en.wikipedia.org/wiki/Canadian_postal_abbreviations_for_provinces_and_territories}{postal abbreviations}
(e.g. "Ontario" would be "ON").

The \code{breaks} argument is slightly buggy at the moment, possibly 
due to how \R tries to split \code{Date} objects. \code{breaks} can 
be either an integer, in which case it will attempt to create that 
many levels (i.e. setting \code{breaks=3} should split the data into 
three date 'blocks'.) \code{breaks} can also be a character vectors, 
such as \code{"quarter"} or \code{"year"} which attempts to split 
the date information accordingly. See \code{\link{cut.Date}} for 
more details and a complete specification of what is allowed for
\code{breaks}.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
The file used to create the map of Canada is from
\href{http://www.geobase.ca/geobase/en/index.html}{GeoBase} and is 
licensed under the \href{http://data.gc.ca/eng/open-government-licence-canada}{Open
Government License - Canada}.
}
\author{
Wen Chen and Joshua Simpson.
}
%\note{

%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, meta)

\dontrun{
group.spatial(ITS1, meta, date.col="Harvestdate", 
              province.col="Province", rank="p", 
              group=c("Ascomycota", "Basidiomycota"), 
              breaks=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
