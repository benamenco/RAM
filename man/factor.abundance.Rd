\name{factor.abundance}
\alias{factor.abundance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Abundance of OTUs by Classification at a Given Taxonomic Rank
For Each Level of A Metadata Category Variable.
}
\description{
This function consumes an OTU, and a rank, as well as various optional
parameters. It creates a stacked bar plot showing the abundance of all
classifications at the given taxonomic rank for each level of a 
metadata category variable.
}
\usage{
factor.abundance(data, rank, top=NULL, count=FALSE, meta=meta, 
                 meta.factor="", drop.unclassified=FALSE, 
                 file=NULL, ext=NULL, height=8, width=16, main="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables with names.
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
  \item{meta}{
metadata.
}
  \item{meta.factor}{
a category variable in metadata
}
  \item{drop.unclassified}{
logical. Should unclassified samples be excluded from the data?
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
  \item{main}{
the title of the plot
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
Wen Chen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(ITS1, ITS2, meta)
data=list(ITS1=ITS1, ITS2=ITS2)

# plot the relative abundance at the class level to the screen, 
# ignoring the unclassified
factor.abundance(data=data, rank="family", meta=meta, 
                 meta.factor=c("Crop"), top=20, 
                 drop.unclassified=TRUE)
\dontrun{
# plot the count abundance at the phylum level to path.tiff
factor.abundance(data=data, rank="family", meta=meta, 
                 meta.factor=c("Crop"), top=20, count=FALSE, 
                 drop.unclassified=TRUE, main="", 
                 file="path/to/tiff", ext="tiff", height=8, width=12)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
