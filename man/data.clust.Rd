\name{data.clust}
\alias{data.clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Hierarchical Cluster Of Samples Based on OTU Table or Taxonomic
Abundance Matrix
}
\description{
This function plot hierarchical cluster Of ecology data set.
}
\usage{
data.clust(data, is.OTU=TRUE, meta, rank=NULL, top=NULL, 
           mode="number", group=4, data.trans=NULL, 
           dist=NULL, clust=NULL, type=NULL, main=NULL, 
           file=NULL, ext=NULL, width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an ecology data set to be analyzed.
}
  \item{is.OTU}{
logical. If an OTU table was provided, \code{is.OTU} should be set 
as \code{TRUE}; otherwise, it should be set as \code{FALSE}.
}
  \item{meta}{
the metadata table associated with ecology data set.
}
  \item{rank}{
optional. If no \code{rank} was provided, the data will be used as
it is, if \code{rank} is provided, if data is an OTU table, it will
be converted to taxonomic abundance matrix at the given rank, no 
change will be made for a data that has already been a taxonomic 
abundance matrix. 
See also \code{\link{tax.abund}} and \code{\link{data.revamp}}
}
  \item{top}{
the top otuIDs or taxa to be considered for the clustering analysis.
See also \code{\link{data.revamp}}
}
  \item{mode}{
either be "number" or "percent". See also \code{\link{data.revamp}}
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
the given method; otherwise use vegdist default Bray-Curtis method. 
See also \code{\link{vegdist}} and \code{\link{gowdis}}.
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
This function returns a tree plot of the hierarchical cluster of 
the samples based on ecological data. 
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

%\seealso{
% ~~further links~~
%}
\examples{
data(ITS1, meta)
\dontrun{
data.clust(data=ITS1, is.OTU=TRUE, data.trans="total", dist="bray", 
           type="fan", meta=meta, group="Plots")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
