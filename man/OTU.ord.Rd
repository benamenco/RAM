\name{OTU.ord}
\alias{OTU.ord}
\alias{OTU.cca}
\alias{OTU.rda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordination Plot For OTUs Using CCA or RDA Analysis
}
\description{
This function consumes an OTU table, metadata factors, and
graphing options, then produces a plot showing the
\code{\link[vegan]{cca}} or \code{\link[vegan]{rda}} analysis 
of the OTU table.
}
\usage{
OTU.ord(otu, meta=meta, factors=NULL, group=NULL,
        na.action=c("na.fail", "na.omit", "na.exclude"),
        rank="g", taxa=NULL, data.trans="total",
        plot.species=TRUE, plot.scaling=-1,
        biplot.scale=NULL, biplot.sig=NULL, biplot.label= TRUE,
        mode=c("rda", "cca"), choice=c(1,2),
        main="", cex.point=3, cex.leg=12, cex.bp=3,
        file=NULL, ext=NULL, width=10, height=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu}{
the OTU table to be used.
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
  \item{na.action}{
choice of one of the following: "na.fail", "na.omit" or 
"na.exclude", see \code{na.action} in \code{\link[vegan]{cca}} 
for detail.
}
  \item{rank}{
the rank to select the taxon groups at.
}
  \item{taxa}{
an integer or a character vector of taxa names at the given 
rank. If it's an integer, will display the top most abundant 
taxa; if it's a vector of taxa names, will plot the selected 
taxa.
}
\item{data.trans}{
a method used to standardize the OTU table. One of 
"\code{total}", "\code{max}", "\code{freq}", "\code{normalize}",
"\code{range}", "\code{standardize}", "\code{pa}", 
"\code{chi.square}", "\code{hellinger}" or "\code{log}"
(see ?decostand).
}
\item{plot.species}{
whether plot sites or taxa, should be reflex to 
\code{plot.scaling}
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
significance cutoff for biplot to be displayed. Currently 
disabled because in the function, calculated ordination model 
cannot be passed to anova test.
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
\item{cex.leg}{
size of lengend name
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
eigenvalues (i.e. from all axes) is explained by that (single) 
axis.
}
\value{
return a list of following: 1) ggplot object; 2) ordination 
model; 3) commodity data and 4) metadata used for the ordination
model.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
%\note{
%
%}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~
\seealso{
\code{\link{decostand}}, \code{\link{Taxa.ord}},
\code{\link{pcoa.plot}}
}
\examples{
data(ITS1, meta)
its1<- filter.OTU(data=list(ITS1=ITS1), percent=0.001)[[1]]
factors=c("City", "Crop", "Harvestmethod", "Ergosterol_ppm")
\dontrun{
# plot sites
ord1 <- OTU.ord(its1,  meta=meta, data.trans="total",
        factors=factors, mode="cca", biplot.sig=0.1,
        taxa=20, biplot.scale=1.5, cex.point=5,
        plot.species=FALSE, rank="f", plot.scaling=3,
        group=c(City="City", Crop="Crop"))
# plot taxa
ord2 <- OTU.ord(its1,  meta=meta, data.trans="total", 
        plot.scaling=-1,
        factors=factors, mode="cca", biplot.sig=0.1,
        taxa=20, biplot.scale=3, cex.point=3,
        plot.species=TRUE, rank="g")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
