\name{group.abundance}
\alias{group.abundance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Abundance of OTUs by Classification at a Given 
Taxonomic Rank
}
\description{
This function consumes an OTU, and a rank, as well as various 
optional parameters. It creates a stacked bar plot showing the 
abundance of all classifications at the given taxonomic rank 
for each sample.
}
\usage{
group.abundance(data, rank,
                top=NULL, count=FALSE, drop.unclassified=FALSE,
                cex.x=NULL, main=NULL, file=NULL, ext=NULL,
                height=8, width=16, bw=FALSE, ggplot2=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables.
}
  \item{rank}{
the taxonomic rank to use.
See \code{\link{RAM.rank.formatting}}.
}
  \item{top}{
the number of groups to select, starting with the most abundant.
If \code{NULL}, all are selected.
}
  \item{count}{
logical. If \code{TRUE}, the numerical counts for each OTU will 
be shown; otherwise the relative abundance will be shown.
}
  \item{drop.unclassified}{
logical. Should unclassified samples be excluded from the data?
}
  \item{cex.x}{
optional. The size of x axis names.
}
  \item{main}{
the title of the plot
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
  \item{bw}{
logical. Should the image be created in black and white?
}
  \item{ggplot2}{
logical. Should the \code{ggplot2} package be used to produce 
the plot, or should the base graphics be used? 
(see ?RAM.plotting).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
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
%% ~Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, ITS2)
# plot the relative abundance at the class level to the screen, 
# ignoring the unclassified
group.abundance(data=list(ITS1=ITS1), rank="phylum",
                drop.unclassified=TRUE)
\dontrun{
# plot the count abundance at the phylum level to path.tiff
group.abundance(data=list(ITS1=ITS1, ITS2=ITS2), rank="g",
                top=10, count=FALSE, drop.unclassified=TRUE,
                main="", file=NULL, ext=NULL,
                height=8, width=16, bw=FALSE, ggplot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
