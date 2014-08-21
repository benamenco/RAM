\name{group.abundance}
\alias{group.abundance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Abundance of OTUs by Classification at a Given Taxonomic Rank
}
\description{
This function consumes an OTU, and a rank, as well as various optional
parameters. It creates a stacked bar plot showing the abundance of all
classifications at the given taxonomic rank for each sample.
}
\usage{
group.abundance(otu1, otu2=NULL, rank, 
                top=NULL, count=FALSE, drop.unclassified=FALSE,
                file=NULL, ext=NULL, labels=c("ITS1", "ITS2"),
                height=8, width=16, bw=FALSE, ggplot2=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be used.
}
\item{otu2}{
the second OTU table to be used.
}
  \item{rank}{
the taxonomic rank to use (see ?RAM.rank.formatting for formatting details).
}
  \item{top}{
the number of groups to select, starting with the most abundant. If \code{NULL},
all are selected.
}
  \item{count}{
logical. If \code{TRUE}, the numerical counts for each OTU will be shown;
otherwise the relative abundance will be shown.
}
  \item{drop.unclassified}{
logical. Should unclassified samples be excluded from the data?
}
  \item{file}{
the file path where the image should be created (see ?RAM.plotting).
}
  \item{ext}{
the file type to be used; one of \code{"pdf"}, \code{"png"}, \code{"tiff"}, 
\code{"bmp"}, \code{"jpg"}, or \code{"svg"}.
}
\item{labels}{
a character vector giving the labels for the panels of the plot.
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
logical. Should the \code{ggplot2} package be used to produce the plot, or
should the base graphics be used? (see ?RAM.plotting).
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, ITS2)

# plot the relative abundance at the class level to the screen, ignoring the
# unclassified
group.abundance(ITS1, rank="class", drop.unclassified=TRUE)

\dontrun{
# plot the count abundance at the phylum level to path.tiff
group.abundance(ITS1, ITS2, rank="p", file="my/file/path", ext="tiff",
		count=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
