\name{phylo_taxonomy}
\alias{phylo_taxonomy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Hierarchical Taxonomic Tree with Relative Abundance
}
\description{
This function plots hierarchical taxonomic tree with relative 
abundance of all taxa at a give rank if category variables are 
provided. Nonsplitting nodes are collapsed as ape::plot.phylo does.
}
\usage{
phylo_taxonomy(otu, rank="order", rank.sep="; ", meta, factors, 
               plot.type="phylogram", edge.width=1, cex=0.7, 
               font = 1, x.lim = NULL, tip.offset=0, tip.cex=0.5,
               thermo=FALSE, thermo.horiz=TRUE, thermo.width=0.5,
               thermo.height=1, node.frame="r", node.bg="white",
               node.col="black", node.width=0.5, node.height=0.6,
               node.cex=0.6, node.font=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu}{
an otu table.
}
  \item{rank}{
taxonomic ranks, see ?RAM.rank.formatting
}
  \item{rank.sep}{
the delimiter that separate the taxonomic ranks in the otu table;
default is "; " (semi colon and a white space) 
}
  \item{meta}{
the metadata table associated with otu table, samples must be in 
the same order in both otu table and metadata.
}
  \item{factors}{
the metadata variables, must be categories
}
  \item{plot.type}{
tree type, default is phylogram. see ?plot.phylo
}
   \item{edge.width}{
see ?ape::plot.phylo
}
  \item{cex}{
size of tip labels.
}
  \item{font}{
font of tip labels
}
  \item{x.lim}{
a numeric vector of length one or two giving the limit(s) of
the x-axis. see ?ape::plot.phylo
}
  \item{tip.offset}{
the distance between tips of the phylogeny and their corresponding labels, see ?plot.phylo
}
  \item{tip.cex}{
size of the pies 
}
  \item{thermo}{
add pies or thermometers to the tips 
}
  \item{thermo.horiz}{
orientiation of the thermometers 
}
  \item{thermo.width}{
width of the thermometers  
}
  \item{thermo.height}{
height of the thermometers  
}
  \item{node.frame}{
type of frame around the nodes  
}
  \item{node.bg}{
background color of text frames of nodes   
}
  \item{node.col}{
color of the nodes  
}
  \item{node.width}{
width of the text frames of nodes   
}
  \item{node.height}{
height of the text frames of nodes   
}
  \item{node.cex}{
size of text of nodes
}
  \item{node.font}{
font of text of nodes 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%This function plots hierarchical taxonomic tree with relative 
%%abundance of all taxa at a give rank if category variables are 
%%provided. Nonsplitting nodes are collapsed as ape::plot.phylo does.
%}
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
ITS1.1<-filter.OTU(data=list(ITS1=ITS1), percent=0.01)[1]
factors=c("Crop", "City")
res<-phylo_taxonomy(otu=ITS1, meta=meta, factors=factors,    
                    rank="order", rank.sep="; ", tip.offset=0, 
                    x.lim=NULL, thermo=FALSE, cex=0.5, 
                    tip.cex=0.5)
names(res)
require("plotKML")
par(mfrow=c(1,2))
for (i in 1:length(factors)) {
  display.pal(res[[1]][[i]])
}
par()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
