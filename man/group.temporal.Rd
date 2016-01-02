\name{group.temporal}
\alias{group.temporal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Temporal Trends for Metadata and Taxon Groups
}
\description{
This function consumes an OTU table and its associated metadata,
and creates a plot showing how the collections of taxonomic 
groups, as well as metadata factors, evolve over time.
}
\usage{
group.temporal(data, meta, date.col, factors, rank, group,
               file = NULL, ext = NULL, height = 8, width = 12)
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
a character vector specifying which column of the metadata has
date information (see \code{\link{RAM.dates}}).
}
  \item{factors}{
a named character vector specifying the names of the metadata
columns to be plotted with the taxon group data.
(see \code{\link{RAM.factors}}). NOTE: these factors must be
\emph{numeric} variables.
}
  \item{rank}{
a character vector specifying the rank of the desired taxon 
groups. Note that all groups should come the same rank.
(see \code{\link{RAM.rank.formatting}}).
}
  \item{group}{
a character vector giving the names of the groups to be plotted.
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
The image created will contain several plots. It will always 
contain a large panel showing the counts collected for the 
specified taxon groups over time, and above that panel (on a 
common x-axis) will be a line graph for each metadata factor 
specified.
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
Wen Chen and Joshua Simpson
}
\note{
If your data has collections being taken roughly annually, you
may have a large amount of "empty space" in the middle of your
plot. Consider subsetting the data by year, and plotting each 
year individually using this function.
}
%% ~Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, meta)
group.temporal(ITS1, meta, date.col="Harvestdate",
               factors=c(Ergosterol="Ergosterol_ppm"),
               rank="p", group=c("Ascomycota", "Basidiomycota"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
