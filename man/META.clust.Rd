\name{META.clust}
\alias{META.clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Hierarchical Cluster Of Metadata
}
\description{
This function plot hierarchical cluster Of metadata.
}
\usage{
META.clust(meta, group=4, data.trans=NULL, dist=NULL, clust=NULL, 
          type=NULL, main="", file=NULL, ext=NULL, width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta}{
the metadata table to be clustered.
}
  \item{group}{
an integer or a metadata variable. If an integar, will cut tree into
corresponding groups and color them accordingly; if a metadata 
variable was provided, tree leaves (sampleIDs) will be colored by 
each level.   
}
  \item{data.trans}{
optional. If was provided, numeric data will be transformed. 
See also \code{\link{decostand}}
}
  \item{dist}{
optional. If was provided, distance matrix will be calculated using 
the give method for numeric variables; otherwise use 
\code{\link{vegdist}} default Bray-Curtis method. If metadata 
include qualitative variables, distance matrix will be calculated 
by \code{\link{gowdis}}. 
}
  \item{clust}{
optional. If was not provided, will use the default agglomeration 
method used by hclust, i.e. "complete". Otherwise, will used user
defined method for clustering. See also \code{\link{hclust}}.
}
  \item{type}{
optional. Can be one of the following: "triangle", "rectangle",
"phylogram", "cladogram", "fan", "unrooted", "radial".
}
  \item{main}{
The title of the plot.  
}
  \item{file}{
optional. Filename that the plot to be saved to.  
}
  \item{ext}{
optional. File type that the plot to be saved to.  
}
  \item{width}{
an integer, width of the plot.  
}
  \item{height}{
an integer, height of the plot.   
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function return a plot of the hierarchical cluster analysis 
on a set of metadata.

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
% \note{
%% ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vegdist}} and \code{\link{gowdis}}.
}
\examples{
data(meta)
  
META.clust(meta=meta, type="fan")
META.clust(meta=meta, type="fan", group="City")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
