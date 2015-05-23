\name{phylog_taxonomy}
\alias{phylog_taxonomy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Hierarchical Taxonomic Tree
}
\description{
This function plots hierarchical taxonomic tree, the leaves are 
taxa at a give rank, nonsplitting nodes are not collapsed as 
ape::plot.phylo does.
}
\usage{
phylog_taxonomy(otu, rank="order", rank.sep="; ", meta, 
                factors=NULL, sel.taxon=NULL, sel.rank=NULL, 
                root="root", cleaves=1, cnodes=1, 
                clabel.leaves = 0.5, clabel.nodes = 0.5, 
                f.phylog = 0.5, sub = TRUE, csub = 1.25, 
                possub = "bottomleft", draw.box = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu}{
an otu table.
}
  \item{rank}{
optional. If no \code{rank} was provided, the data will be used as
it is, if \code{rank} is provided, if data is an OTU table, it will
be converted to taxonomic abundance matrix at the given rank, no 
change will be made for a data that has already been a taxonomic 
abundance matrix. 
See also \code{\link{tax.abund}} and \code{\link{data.revamp}}
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
  \item{sel.taxon}{
optional, a selected taxon at higher taxonomic rank (i.e. at 
phylum level). If given, only descendents of this taxa will be plotted.
}
   \item{sel.rank}{
the rank of \code{sel.taxon}
}
  \item{root}{
optional, the name of the root of the tree. If not given, will use "root"
}
  \item{cleaves}{
see ?ade4::plot.phylog
}
  \item{cnodes}{
see ?ade4::plot.phylog
}
  \item{clabel.leaves}{
see ?ade4::plot.phylog
}
  \item{clabel.nodes}{
see ?ade4::plot.phylog 
}
  \item{f.phylog}{
see ?ade4::plot.phylog
}

  \item{sub}{
see ?ade4::plot.phylog
}
  \item{csub}{
see ?ade4::plot.phylog
}
  \item{possub}{
see ?ade4::plot.phylog 
}
  \item{draw.box}{
see ?ade4::plot.phylog
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%This function plots hierarchical taxonomic tree, the leaves are 
%%taxa at a give rank, nonsplitting nodes are not collapsed as 
%%ape::plot.phylo does.
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
ITS1.1<-filter.OTU(data=list(ITS1=ITS1), percent=0.01)[[1]]
factors=c("Crop", "City")
phylog_taxonomy(otu=ITS1.1, rank="family", rank.sep="; ", 
                meta=meta, factors=NULL, sel.taxon=NULL,
                sel.rank=NULL, cleaves=1, cnodes=.5, root="k__Fungi",
                clabel.leaves = 0.5, clabel.nodes = 0.5, 
                f.phylog = 0.8, sub = TRUE, csub = 1.25, 
                possub = "bottomleft", draw.box = TRUE)

phylog_taxonomy(otu=ITS1.1, rank="family", rank.sep="; ", 
                meta=meta, factors=c("Crop", "City", "Plots"),
                sel.taxon=NULL, sel.rank=NULL, cleaves=1, 
                cnodes=.5, clabel.leaves = 0.5, clabel.nodes = 0.5,
                f.phylog = 0.8, sub = TRUE, csub = 1.25, 
                possub = "bottomleft", draw.box = TRUE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
