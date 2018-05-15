\name{network_data}
\alias{network_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates Nodes and Edge-List For An OTU Table
}
\description{
This function creates nodes and edge-list for an otu table, 
which can be used for network plotting
}
\usage{
network_data(data, is.OTU=TRUE, metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an otu table of a species abundance matrix
}
  \item{is.OTU}{
whether or not \code{data} is an otu table
}
  \item{metadata}{
associated metadta of the otu table, should have same sampleIDs
}
}
%\details{
%%
%}
\value{
Return a list of network nodes and edges based on an otu table 
and the associated metadata
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
%\seealso{
% ~~further links~~
%}
\examples{
data(ITS1, meta)
\dontrun{
ITS1.01<-filter.OTU(data=list(ITS1=ITS1), percent=0.01)[[1]]

# create nodes and edges lists
b<-network_data(ITS1.01, is.OTU=TRUE, meta)
b_node<-b[[1]]
b_edge<-b[[2]]
head(b_node)
head(b_edge)

library(igraph)
b1<-graph.data.frame(b_edge, directed=FALSE)
lev <- levels(factor(E(b1)$Crop))

# vertices size
V(b1)$size<-degree(b1)

# vertices color
Crop1<-rownames(meta)[meta$Crop=="Crop1"]
Crop2<-rownames(meta)[meta$Crop=="Crop2"]

## vertices representing samples from crop1 will be in red,
## vertices representing samples from crop2 will be in blue;
## vertices representing taxa will be in pink
V(b1)$color[which(is.element(V(b1)$name, Crop1))]<-"red"
V(b1)$color[which(is.element(V(b1)$name, Crop2))]<-"blue"
V(b1)$color[-which(is.element(V(b1)$name, c(Crop1, Crop2)))]<-"pink"
V(b1)$color

# edges color
col<-c("red", "blue")
for (i in 1:length(lev) ) { 
  E(b1)$color[E(b1)$Crop==lev[i]] <- col[i]
}

# plot
plot(b1, vertex.label.font=2,
         vertex.label.cex=0.5,
         layout=layout.kamada.kawai)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
