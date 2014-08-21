\name{group.heatmap}
\alias{group.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot OTU Abundance at a Given Rank with Metadata Annotation
}
\description{
This function plots the abundance for all taxon groups at a given rank.
Additionally, it can display metadata for the samples as annotations along the
rows of the heatmap.
}
\usage{
group.heatmap(data, meta, rank, factors, 
              top=NULL, count=FALSE, drop.unclassified=FALSE,
              cut=NULL, file=NULL, ext=NULL, width=9, height=9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be used.
}
  \item{meta}{
the metadata table to be used.
}
  \item{rank}{
the taxonomic rank to use (see \code{\link{RAM.rank.formatting}} for formatting
details).
}
  \item{factors}{
a named character vector indicating the columns of the metadata table to be used 
(see \code{\link{RAM.factors}}).
}
  \item{top}{
the number of groups to select, starting with the most abundant. If \code{NULL},
all are selected.
}
\item{count}{
logical. If \code{TRUE}, the numerical counts for each group will be shown;
otherwise the relative abundance will be shown.
}
  \item{drop.unclassified}{
logical. Should OTUs labelled "unclassified" or missing classification at the
given taxonomic rank be excluded?
}
  \item{cut}{
the height at which to cut the sample tree, this will create distinct coloured 
groups. Currently this will allow for at most nine groups (see Details).
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
A warning from \code{brewer.pal} indicating "n too large, allowed maximum for 
palette Pastel1 is 9" means that the cut height is too low to allow for that
many groups. This should be fixed in a future release.
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
\code{\link{decostand}}
}
\examples{
data(ITS1, meta)

group.heatmap(ITS1, meta, rank="c", factors=c(City="City", Plot="Plots"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
