\name{Taxa.ord}
\alias{Taxa.ord}
\alias{Taxa.cca}
\alias{Taxa.rda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordination Plot For Taxa Groups Using CCA or RDA Analysis
}
\description{
This function consumes an ecology data set, metadata factors, and 
graphing options, then produces a plot showing the vegan::cca or
vegan::rda analysis.
}
\usage{
Taxa.ord(data, is.OTU=TRUE, meta=meta, factors=NULL, 
         group=NULL, rank="g", taxa=10, data.trans="total", 
         plot.species=TRUE, plot.scaling=-1, 
         biplot.scale=NULL, biplot.sig=NULL, biplot.label= TRUE, 
         mode=c("rda", "cca"), choice=c(1,2), main="", 
         cex.point=3, cex.label=1, cex.leg=12, cex.bp=3, cex.text=3,
         file=NULL, ext=NULL, width=10, height=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an ecology data set, either an otu table or a taxonomy abundance 
matrix.
}
  \item{is.OTU}{
whether or not the data an otu table
}
  \item{meta}{
the metadata table to be used.
}
  \item{factors}{
a named character vector of length 1 or 2 specifying metadata 
factors for the samples in the OTU table (see Details).
}
  \item{group}{
a named character vector of length 1 or 2 specifying metadata 
factors for the samples in the OTU table (see Details).
}
  \item{rank}{
the rank to select the taxon groups at. 
}
  \item{taxa}{
an integer or a character vector of taxa names at the given rank.
if \code{integer}, plot the top most abundant taxa, otherwise plot 
the taxa in the vector. 
}
  \item{data.trans}{
a method used to standardize the OTU table. One of "\code{total}", 
"\code{max}", "\code{freq}", "\code{normalize}", "\code{range}", 
"\code{standardize}", "\code{pa}", "\code{chi.square}", 
"\code{hellinger}" or "\code{log}" 
(see ?decostand).
}
  \item{plot.species}{
whether plot sites or taxa, should be reflex to \code{plot.scaling}
}
  \item{plot.scaling}{
one of the following: 1, 2, 3, or -1. 
See \code{scaling} in \code{\link{plot.cca}} for detail.
See also \code{\link{ordiplot}}
}
  \item{biplot.scale}{
a numeric number, length of the biplot arrows
} 
\item{biplot.sig}{
significance cutoff for biplot to be displayed. Currently disabled 
because in the function, ordination model calcuated cannot be passed
to anova test.
}
  \item{biplot.label}{
whether or not to plot biplot
} 
  \item{mode}{
one of the following: "cca" or "rda".
} 
\item{choice}{
the chosen axes
} 
  \item{main}{
title of the plot
} 
  \item{cex.point}{
size of points
} 
  \item{cex.label}{
size of taxa lables
}
  \item{cex.leg}{
size of lengend name
} 
  \item{cex.text}{
size of taxon names if \code{plot.species} is set \code{TRUE}
}
  \item{cex.bp}{
size of biplot labels
}
  \item{file}{
the file path where the image should be created 
(see ?RAM.plotting).
}
  \item{ext}{
the file type to be used; one of \code{"pdf"}, \code{"png"}, 
\code{"tiff"}, \code{"bmp"}, \code{"jpg"}, or \code{"svg"}.
}
  \item{width}{
the width of the image to be created (in inches).
}
  \item{height}{
the height of the image to be created (in inches).
}
}
\details{
\code{group} should be a named character vector specifying the 
names of the columns to be used from \code{meta} 
(see \code{\link{RAM.factors}}). 
The values on the axes denote what fraction of the sum of all 
eigenvalues (i.e. from all axes) is explained by that (single) axis.
}
\value{
return a list of following: 1) ggplot object; 2) ordination model; 
3) commodity data and 4) metadata used for the ordination model.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
\note{
The labels for the taxa points are placed above, below, or next 
to the point itself at random. If labels are outside of the 
plotting area, or overlapping with each other, run your command 
again (without changing any arguments!) and the labels should move 
to new positions. Repeat until they are placed appropriately. This 
is done to ensure even tightly-grouped samples, or samples near 
the edge of the plot, have their labels shown. If the labels are 
too distracting, remember that they can be turned off by setting 
\code{plot.species = FALSE}.
}

%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link{decostand}}, \code{\link{OTU.ord}}, 
\code{\link{pcoa.plot}}
}
\examples{
data(ITS1, meta) 
its1<- filter.OTU(data=list(ITS1=ITS1), percent=0.001)[[1]]
factors=c("City", "Crop", "Harvestmethod", "Ergosterol_ppm")                                        
\dontrun{
ord <- Taxa.ord(its1,  meta=meta, data.trans="total", 
        factors=factors, mode="cca", biplot.sig=0.1, 
        taxa=20, biplot.scale=1.5, cex.point=5, cex.label=1,
        plot.species=TRUE, rank="g", plot.scaling=3,
        group=c(City="City", Crop="Crop"), biplot.label=FALSE)
names(ord)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
