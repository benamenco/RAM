\name{top.groups.plot}
\alias{top.groups.plot}
\alias{group.top.percent}
\alias{group.top.number}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Top Taxon Groups
}
\description{
These functions consume two OTU tables, along with (optionally) a 
file name and a parameter \code{top}. They create a box plot of 
the \code{top} number of OTUs (for \code{plot.top.number}), or all 
OTUs with relative abundance above \code{top} percent 
(for \code{plot.top.percent}) at the taxonomic ranks phylum,
class, order, family, and genus.
}
\usage{
group.top.number(data, top=10, ranks=c("p","c","o","f","g"), 
                 drop.unclassified=FALSE, cex.x=NULL,
                 main=NULL,file=NULL, ext=NULL, height=8, width=16, 
                 bw=FALSE, ggplot2=TRUE)
                
group.top.percent(data, top=10, ranks=c("p","c","o","f","g"),
                  drop.unclassified=FALSE, cex.x=NULL,
                  main=NULL, file=NULL, ext=NULL, height=8, 
                  width=16, bw=FALSE, ggplot2=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables.
}
\item{top}{
the number of OTUs to select (for \code{top.number}), or the minimum 
relative abundance threshold to use for selecting OTUs 
(for \code{top.percent}).
}
  \item{ranks}{
a vector of the taxonomic ranks.
See also \code{\link{RAM.rank.formatting}}
}
  \item{drop.unclassified}{
logical. Should OTUs labelled "unclassified" or missing 
classification at the given taxonomic rank be excluded?
}
\item{cex.x}{
optional. The size of the x axis names.
}
\item{main}{
the title of the plot
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
  \item{bw}{
logical. Should the image be created in black and white?
}
  \item{ggplot2}{
logical. Should the \code{ggplot2} package be used to produce the 
plot, or should the base graphics be used? (see ?RAM.plotting).  
}
}

\note{
Please be aware that the 'whiskers' in the plot may differ depending 
on the setting of \code{ggplot2}. 
Please see \code{\link{geom_boxplot}} \code{\link{boxplot}}, and
\code{\link{boxplot.stats}} for more information on how the 
whiskers are calculated.
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Wen Chen and Joshua Simpson.
}

\seealso{
\code{\link{RAM.plotting}}
}
\examples{
\dontrun{
data(ITS1, ITS2)

# plots the top 10 OTUs (by default) at five ranks
group.top.percent(data=list(ITS1=ITS1, ITS2=ITS2), top=10) 

# plots all OTUs w/ relative abundance > 10\% (as specified) at 
# five ranks and saves the result as a .tiff file 
# (only using ITS1 data)
group.top.percent(data=list(ITS1=ITS1, ITS2=ITS2), top=10, 
                 file="my/file/path", ext="tiff")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
