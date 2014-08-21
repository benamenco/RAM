\name{dissim.plot}
\alias{dissim.plot}
\alias{dissim.clust.plot}
\alias{dissim.eig.plot}
\alias{dissim.alleig.plot}
\alias{dissim.ord.plot}
\alias{dissim.GOF.plot}
\alias{dissim.tree.plot}
\alias{dissim.pvar.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Dissimilarity Matrix Data for Different Methods
}
\description{
These functions all produce a plot of some measure related to dissimilarity
matrices. All of these functions allow you to specify a vector of methods to be
used when creating the plot.
}
\usage{
dissim.clust.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL,
                  clust.methods=NULL)
dissim.eig.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL)
dissim.alleig.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL)
dissim.ord.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL, k=NULL)
dissim.GOF.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL)
dissim.tree.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL,
                 clust.methods=NULL)
dissim.pvar.plot(otu1, otu2=NULL, file=NULL, dist.methods=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the otu1 data to be used.
}
  \item{otu2}{
the otu2 data to be used.
}
  \item{file}{
the file path for the plot. If not provided (defaults to NULL), then the plot is
displayed to the screen. If \code{file} is provided, that is where the .tiff
file will be created.
}
  \item{dist.methods}{
a character vector representing the dissimilarity indices to be used; each 
element must be one of one of "\code{manhattan}", "\code{euclidean}", 
"\code{canberra}", "\code{bray}", "\code{kulczynski}", "\code{jaccard}", 
"\code{gower}", "\code{altGower}","\code{morisita}", "\code{horn}", 
"\code{mountford}", "\code{raup}", "\code{binomial}","\code{chao}", or 
"\code{cao}".   
}
  \item{clust.methods}{
a character vector representing the methods used for clustering the data. 
Each element must be one of "\code{ward}", "\code{single}", "\code{complete}",
"\code{average}", "\code{mcquitty}", "\code{median}", or "\code{centroid}".   
}
  \item{k}{the number of dimensions desired. If \code{NULL}, the maximum value will be 
calculated and used.}
}
\details{
All of these functions (other than \code{dissim.alleig.plot}) call their
\code{dissim.X} counterparts and plot the data. If \code{file} is given, a
\code{.tiff} file will be created at \code{file}; otherwise the plot is
displayed to the screen.
}
\value{
All functions create a plot and return the plotted data invisibly.
\item{\code{dissim.clust.plot }}{plots a hierarchical clustering of the
dissimilarity matrix.}
\item{\code{dissim.eig.plot }}{plots a bar plot of the eigenvalues of the
dissimilarity matrix.}
\item{\code{dissim.alleig.plot }}{plots a line plot showing the relative
importance of all eigenvalues for a variety of methods.}
\item{\code{dissim.ord.plot }}{plots a scatter plot comparing the "euclidean"
distances among all samples in ordination space to the dissimilarity matrix distances.}
\item{\code{dissim.GOF.plot }}{plots a scatter plot of the goodness of fit
values of the dissimilarity matrix, for various numbers of dimensions used.}
\item{\code{dissim.tree.plot }}{plots a scatter plot comparing the tree
distances to the dissimilarity matrix distances.}
\item{\code{dissim.pvar.plot }}{plots a bar plot showing the percent variation
explained by each axis (where each sample corresponds to an axis).}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson
}
\note{
If \code{file} does not end in ".tiff", then ".tiff" will be appended to the end
of \code{file}.

Function \code{dissim.alleig.plot} uses the \code{ggplot2} package for creating
the plot, and returns the plot object. This means that you can store this plot
and add other features manually, if desired (see Examples).
}

\seealso{
\code{\link{vegdist}}, \code{\link{hclust}}, \code{\link{dissim}},
\code{\link{ggplot}}
}
\examples{
data(ITS1, ITS2)

# show percent variation for only ITS1 with default methods
dissim.pvar.plot(ITS1)

# show clustering for ITS1 and ITS2 for set methods
dissim.clust.plot(ITS1, ITS2, dist.methods=c("morisita", "bray"), 
                  clust.methods=c("average", "centroid"))

\dontrun{
# dissim.alleig.plot returns a ggplot2 object:
my.eig.plot <- dissim.alleig.plot(ITS1, ITS2)
class(my.eig.plot) # returns "gg" "ggplot"
my.eig.plot # view the plot
# update the title, then view the updated plot
my.eig.plot <- my.eig.plot + ggtitle("My New Title")
my.eig.plot

# save an image (named file.tiff) with GOF values for ITS1 and ITS2, using
# default methods
dissim.GOF.plot(ITS1, ITS2, file="~/Documents/my/file")}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
