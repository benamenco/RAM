\name{group.Taxa.bar}
\alias{group.Taxa.bar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Barplot Of Taxa Distribution In Groups
}
\description{
This function do a barplot to show the distribution of selected
taxa in each level of a given metadata variable
}
\usage{
group.Taxa.bar(data, is.OTU=TRUE, rank="g", taxa="", 
               meta, meta.factor="", func="sum", cex.y=5, cex.x=5, 
               cex.main=10, bar.width=NULL, RAM.theme=NULL,
               col.pal=NULL, main="", file=NULL, ext=NULL, 
               height=8, width=16)
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
  \item{meta}{
the metadata table to be used (must have same samples as \code{data}.
}
  \item{meta.factor}{
a character string. Must be one of the metadata variables.
}
  \item{func}{
function to be used to aggregate count data, e.g. "sum" or "mean",
default is "sum"
}
  \item{cex.y}{
size of y axis tick labels.
}
  \item{cex.x}{
size of x axis tick labels.
}
  \item{cex.main}{
size of title.
}
  \item{bar.width}{
width of bars
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
\details{
To use customized color palettes, must pass a vector of color names
or hexadecimals. See examples for detail.
}
\value{
This function returns a barplot.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% notes
%}
\examples{
data(ITS1, ITS2, meta)
taxa <- c("Fusarium", "Alternaria", "Cladosporium")
group.Taxa.bar(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE, 
               rank="g", taxa=taxa, meta=meta, meta.factor="City", 
               cex.y=5, cex.x=5, bar.width=1, RAM.theme=RAM.color())
\dontrun{
# change default color schemes
col <- c("dodgerblue1", "darkcyan")
taxa.1 <- c("Fusarium", "Alternaria", "Cladosporium", "Verticillium",
           "Kondoa")
group.Taxa.bar(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE, 
               rank="g", taxa=taxa.1, meta=meta, meta.factor="City", 
               cex.y=5, cex.x=5, bar.width=1, RAM.theme=NULL,
               col.pal=col)
group.Taxa.bar(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE, 
               rank="g", taxa=taxa.1, meta=meta, meta.factor="City", 
               cex.y=5, cex.x=5, bar.width=1, RAM.theme=NULL,
               col.pal=col, func="mean")
# change ggplot theme
group.Taxa.bar(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE, 
               rank="g", taxa=taxa.1, meta=meta, meta.factor="City", 
               cex.y=5, cex.x=5, bar.width=1, 
               col.pal=col, RAM.theme=RAM.border())
# save the plot
group.Taxa.bar(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE, 
               rank="g", taxa=taxa.1, meta=meta, meta.factor="City", 
               cex.y=5, cex.x=5, bar.width=1, RAM.theme=NULL,
               col.pal=col,main="", file="path/to/filename.pdf", 
               ext="pdf", height=8, width=16)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
