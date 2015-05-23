\name{group.heatmap}
\alias{group.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot OTU Abundance at a Given Rank with Metadata Annotation
}
\description{
This function plots the abundance for all taxon groups at a given 
rank. Additionally, it can display metadata for the samples as 
annotations along the rows of the heatmap.
}
\usage{
group.heatmap(data, is.OTU=TRUE, meta, rank, factors, 
              top=25, remove.unclassified=TRUE, 
              stand.method=NULL,
              dist.method="bray",
              hclust.method="average",
              dendro.row.status="yes",
              dendro.col.status="hidden",
              row.labels=TRUE, row.cex=1,
              cut=NULL, file=NULL, ext=NULL, 
              width=9, height=9)
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
the taxonomic rank to use (see \code{\link{RAM.rank.formatting}} 
for formatting details).
}
  \item{factors}{
a named character vector indicating the columns of the metadata 
table to be used (see \code{\link{RAM.factors}}).
}
  \item{top}{
the number of groups to select, starting with the most abundant. 
If \code{NULL}, all are selected.
}
  \item{remove.unclassified}{
logical. Define whether OTUs labelled "unclassified" or missing 
classification at the given taxonomic rank should be excluded.
}
  \item{stand.method}{
optional. Data standardization methods specified in \code{\link[vegan]{decostand}}.
}
  \item{dist.method}{
distance matrix calculation methods specified \code{\link[vegan]{vegdist}}.
}
  \item{hclust.method}{
the agglomeration methods specified in \code{\link{hclust}}. 
This should be unambiguous abbreviation of one of the following: 
'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 
'median' or 'centroid'.
}
  \item{dendro.row.status}{
whether or not to use the dendrogram to re-order the rows. 
It should be one of the following: 'yes' that use the dendrogram 
to re-order the rows/columns and display the dendrogram; 'hidden'
means re-rorder, but do not display; 'no' means do not use the 
dendrogram at all.
See also \code{\link[Heatplus]{annHeatmap2}}
}
  \item{dendro.col.status}{
whether or not to use the dendrogram to re-order the columns. 
It should be one of the following: 'yes' that use the dendrogram 
to re-order the rows/columns and display the dendrogram; 'hidden'
means re-rorder, but do not display; 'no' means do not use the 
dendrogram at all.
See also \code{\link[Heatplus]{annHeatmap2}}
}
  \item{row.labels}{
logical. Whether or not to plot the row labels.
}
  \item{row.cex}{
the size of row labels if \code{row.labels} is \code{TRUE}
}
  \item{cut}{
the height at which to cut the sample tree, this will create 
distinct coloured groups. Currently this will allow for at most 
nine groups (see Details).
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
\details{
A warning from \code{brewer.pal} indicating "n too large, allowed 
maximum for palette Pastel1 is 9" means that the cut height is 
too low to allow for that many groups. This should be fixed in a 
future release.
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

%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link{decostand}}, \code{\link[Heatplus]{annHeatmap2}}
}

\examples{
data(ITS1, meta)
\dontrun{
library("Heatplus")
library("RColorBrewer")

group.heatmap(ITS1, is.OTU=TRUE, meta=meta, rank="c", factors=c(Crop="Crop", 
              City="City"), stand.method="chi", dist.method="euc", 
              hclust.method="mcquitty", cut=0.5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
