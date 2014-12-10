\name{pcoa.plot}
\alias{pcoa.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a PCoA plot for an OTU Table
}
\description{
This function consumes an OTU table, metadata factors, and 
graphing options, then produces a plot showing the PCoA analysis 
of the OTU table.
}
\usage{
pcoa.plot(data, is.OTU=TRUE, meta, factors, rank, stand.method = NULL,
          dist.method = "morisita", sample.labels = TRUE, top = 20, 
          ellipse = FALSE, main = NULL, file = NULL, ext = NULL, 
          height = 8, width = 10, ggplot2 = TRUE, bw = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an OTU table or taxonomic abundance matrix to be used.
}
  \item{is.OTU}{
logical. Whether or not the input data an OTU table.
}
  \item{meta}{
the metadata table to be used.
}
  \item{factors}{
a named character vector of length 1 or 2 specifying metadata 
factors for the samples in the OTU table (see Details).
}
  \item{rank}{
the rank to select the taxon groups at. 
For an OTU table, if \code{rank} is set \code{NULL}, distance 
matrix will be calcuated using all OTUs, otherwise, the OTU table 
will be transformed to taxonomic abundance matrix before the 
calculation of the distance matrix. If a taxonomic abundance matrix
is provided, i.e. \code{is.OTU} is set \code{TRUE}, then the \code{rank}
will be ignored.
}
\item{stand.method}{
a method used to standardize the OTU table. One of "\code{total}", 
"\code{max}", "\code{freq}", "\code{normalize}", "\code{range}", 
"\code{standardize}", "\code{pa}", "\code{chi.square}", 
"\code{hellinger}" or "\code{log}" 
(see ?decostand).
}
  \item{dist.method}{
the dissimilarity index to be used; one of "\code{manhattan}",
"\code{euclidean}", "\code{canberra}", "\code{bray}", 
"\code{kulczynski}", "\code{jaccard}", "\code{gower}", 
"\code{altGower}","\code{morisita}", "\code{horn}", 
"\code{mountford}", "\code{raup}", "\code{binomial}",
"\code{chao}", or "\code{cao}" 
(see \code{\link{vegdist}}).
}
\item{sample.labels}{
logical. Should the labels for the samples be displayed?
}
  \item{top}{
how many taxon groups should be displayed, starting from the most 
abundant.
}
  \item{ellipse}{
which of the metadata factors (if any) should have ellipses 
plotted around them. Must be one of 1, 2, or FALSE 
(see Details).
}
  \item{main}{
The title of the plot.
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
  \item{ggplot2}{
logical. Should the \code{ggplot2} package be used to produce the 
plot, or should the base graphics be used? (see ?RAM.plotting).  
}
  \item{bw}{
logical. Should the image be created in black and white?
}
}
\details{
This function uses \code{\link{pco}} in the labdsv package for the 
Principal coordinates analysis (PCoA). The distance matrix was 
square rooted before being passed to \code{\link{pco}} to avoid 
negative eigenvalues.
\code{factors} should be a named character vector specifying the 
names of the columns to be used from \code{meta} 
(see \code{\link{RAM.factors}}). 
Those columns should be factors; if they are not, a warning is 
raised and they are coerced to factors (see \code{\link{factor}}). 
A warning is also raised when a factor has more than 9 levels, 
as that is the most colours the current palettes support.

The values on the axes denote what fraction of the sum of all 
eigenvalues (i.e. from all axes) is explained by that (single) axis.

When \code{ellipse = FALSE}, no ellipses will be plotted. When 
\code{ellipse} is a number, that 'number' metadata factor will 
have ellipses plotted. 
For example, if \code{factors = c(Crop="Crop", City="City")} and 
\code{ellipse = 1}, ellipses will be plotted for the different 
crops, but NOT the cities. Setting \code{factors = c(City="City")} 
and \code{ellipse = 2} is invalid, since there is no second 
metadata factor given. Ellipses can only be plotted for one factor 
currently. Furthermore, there need to be at least 3 samples for 
every level in every item in \code{factors}, otherwise ellipses 
cannot be plotted.
}
\value{
When \code{ggplot2 = TRUE}, a ggplot object is returned; otherwise 
nothing ism returned (but the plot is shown on screen).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
\note{
The labels for the sample points are placed above, below, or next 
to the point itself at random. If labels are outside of the 
plotting area, or overlapping with each other, run your command 
again (without changing any arguments!) and the labels should move 
to new positions. Repeat until they are placed appropriately. This 
is done to ensure even tightly-grouped samples, or samples near 
the edge of the plot, have their labels shown. If the labels are 
too distracting, remember that they can be turned off by setting 
\code{sample.labels = FALSE}.
}

%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link{vegdist}}
}
\examples{
data(ITS1, meta)

# The argument for factors is a vector of length two; the first 
# item is # Crop, which is a column from meta, and the second item 
# is City, another # column from meta.

pcoa.plot(ITS1, meta=meta, rank="c", 
          factors=c(Crop="Crop", City="City"))
                              
\dontrun{
# If you want to customize legend labels and plot the top 20 taxon 
# groups at genus:
pcoa.plot(ITS1, meta=meta, rank="g", main="PCoA plot",
          factors=c(Place="City", Harvest_Method="Harvestmethod"))
                                           
# In black & white, using base graphics:
pcoa.plot(ITS1, meta=meta, rank="c", factors=c(Plot="Plots"), 
          ggplot=F, bw=T)

pcoa.plot(ITS1, meta=meta, rank="c", factors=c(Plot="Plots"), 
          ggplot=F, bw=T, dist.method="euc", stand.method="hell")

# Focus on the samples: hide all groups and plot ellipses for Crop:
pcoa.plot(ITS1, meta=meta, rank="g", 
          factors=c(Crop="Crop", City="City"), 
          ellipse=1, sample.labels=FALSE, top=0)

# Standardize the data before calculating distances:
pcoa.plot(ITS1, meta=meta, rank="g", factors=c(City="City"),  
          stand.method="chi.square",
          dist.method="euclidean")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
