\name{group.abund.Taxa}
\alias{group.abund.Taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barplot Of Distribution Of Taxa In Groups
}
\description{
This function do a barplot to show the distribution of selected
taxa in each level of a given metadata variable
}
\usage{
group.abund.Taxa(data, is.OTU=TRUE, rank="g", taxa, 
                 drop.unclassified=FALSE, bar.width=NULL,
                 meta, meta.factor="", RAM.theme=NULL, 
                 col.pal=NULL, main="",
                 file=NULL, ext=NULL, height=8, width=16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of otu tables or taxonomic abundance matrices.
See also \code{\link{RAM.input.formatting}}.
}
  \item{is.OTU}{
logical. If an OTU table was provided, \code{is.OTU} should be set 
as \code{TRUE}; otherwise, it should be set as \code{FALSE}.
}
  \item{rank}{
a single taxonomic rank. See also \code{\link{RAM.rank.formatting}}
}
  \item{taxa}{
a vector containing taxa names for plotting.
}
  \item{drop.unclassified}{
logical. Whether or not drop the unclassified taxon groups.
}
  \item{bar.width}{
width of bars
}
  \item{meta}{
the metadata table to be used (must have same samples as \code{data}.
}
  \item{meta.factor}{
a character string. Must be one of the metadata variables.
}
  \item{RAM.theme}{
customized ggplot_theme in RAM. 
See also ?theme_ggplot.
}
  \item{col.pal}{
color palettes to be used. 
}
  \item{main}{
a character string. The title of the plot, default is an empty string.
}
  \item{file}{
filename to save the plot.
}
  \item{ext}{
filename extension, the type of image to be saved to.
}
  \item{width}{
an integer, width of the plot.
}
  \item{height}{
an integer, height of the plot.
}
}
%\details{
%%  Maybe explain full/restricted models
%}
\value{
This function returns a Barplot of the distribution of seleted 
taxa within each level of a given metadata variable.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
\note{
This funtion provide an alternative view of taxa distribution as
\code{\link{group.Taxa.bar}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% notes
%}
\examples{
data(ITS1, ITS2, meta)
taxa <- c("Fusarium", "Alternaria", "Cladosporium")
group.abund.Taxa(data=list(ITS1=ITS1, ITS2=ITS2), taxa=taxa, 
                 meta=meta, meta.factor="Crop", 
                 drop.unclassified=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
