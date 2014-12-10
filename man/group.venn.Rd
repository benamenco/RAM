\name{group.venn}
\alias{group.venn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Venn Diagram For Two To Five Sets With Item Labels
}
\description{
This function use \code{\link{draw.pairwise.venn}} to creates a 
venn diagram for two vectors
}
\usage{
group.venn(vectors, cat.cex=1.5, cex=1,
           cat.pos=NULL, cat.dist=NULL, 
           label=TRUE, lab.cex=1, 
           lab.col= "black", fill=NULL, 
           file=NULL, ext=NULL, width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vectors}{
a list of vectors with names. See also \code{\link{RAM.input.formatting}}.
}
  \item{cat.cex}{
size of the category names. 
(see \code{\link{venn.diagram}} for details).
}
  \item{cex}{
size of the label of the circles. 
(see \code{\link{venn.diagram}} for details).
}
  \item{cat.pos}{
optional. Location of the category names along the circles.
(see \code{\link{venn.diagram}} for details).
}
  \item{cat.dist}{
optional. Distance of the category names to the circles.
(see \code{\link{venn.diagram}} for details).
}
  \item{label}{
logical. If \code{TRUE}, will plot the item labels for 2 data sets.
For more than 2 datasets or this is set as \code{FALSE}, the labels 
will be numbers for each circle.  
(see \code{\link{venn.diagram}} for details).
}
  \item{lab.cex}{
size of the labels. 
}
  \item{lab.col}{
color of the labels. 
}
  \item{fill}{
optional, color of the circles. 
(see \code{\link{venn.diagram}} for details).
}
  \item{file}{
the file path where the image should be created (see ?RAM.plotting).
}
  \item{ext}{
filename extension, the type of image to be saved to.
}
  \item{height}{
the height of the image to be created (in inches).
}
  \item{width}{
the width of the image to be created (in inches).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{group.venn} returns a venn diagram for 2 to 5 sets. The user 
can choose to place item labels for 2 sets of data, however, the 
label locations can be wrong if the the smaller data set is part of
the bigger data set, in this case, set \code{label} as \code{FALSE}.
If the input datasets is more than 2, \code{label} will be ignored.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
%\note{
%
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
see \code{\link{venn.diagram}}
}
\examples{
data(ITS1, meta)
# core OTUs
core <- core.OTU.rank(data=list(ITS1=ITS1), meta=meta, rank="g", 
                     meta.factor="Crop", percent=1)

# taxa that core OTUs assigned to    
core.Crop1 <- core$ITS1$Crop1$taxa
core.Crop2 <- core$ITS1$Crop2$taxa

# venn plot
vectors <- list(Core_Crop1=core.Crop1, Core_Crop2=core.Crop2)
group.venn(vectors=vectors, label=TRUE, cat.pos=c(330, 150), 
           lab.cex=0.7)
\dontrun{
group.venn(vectors=vectors, label=FALSE, cat.pos=c(330, 150), 
           lab.cex=0.7, cex=3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
