\name{dissim}
\alias{dissim}
\alias{dissim.clust}
\alias{dissim.eig}
\alias{dissim.ord}
\alias{dissim.GOF}
\alias{dissim.tree}
\alias{dissim.pvar}
\title{
Calculate Dissimilarity Matrix Data
}
\description{
These functions calculate different measures related to dissimilarity 
matrices. All of these functions allow you to specify one of many 
dissimilarity indices to be used.
}
\usage{
dissim.clust(data, dist.method="morisita", clust.method="average")
dissim.eig(data, method="morisita")
dissim.ord(data, dist.method="morisita", k=NULL)
dissim.GOF(data, method="morisita")
dissim.tree(data, dist.method="morisita", clust.method="average")
dissim.pvar(data, method="morisita")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be used.
}
  \item{method}{
the dissimilarity index to be used; one of "\code{manhattan}",
"\code{euclidean}", "\code{canberra}", "\code{bray}", "\code{kulczynski}", 
"\code{jaccard}", "\code{gower}", "\code{altGower}","\code{morisita}", 
"\code{horn}", "\code{mountford}", "\code{raup}", "\code{binomial}", 
"\code{chao}", or "\code{cao}".   
}
  \item{dist.method}{
same as \code{method}, but named differently for functions with other 
\code{method} parameters.
}
  \item{k}{
the number of dimensions desired. If \code{NULL}, the maximum value will be 
calculated and used.
}
  \item{clust.method}{
the method used for clustering the data. Must be one of "\code{ward}",
"\code{single}", "\code{complete}", "\code{average}", "\code{mcquitty}", 
"\code{median}", or "\code{centroid}".   
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{\code{dissim.clust }}{returns a hierarchical clustering of the
dissimilarity matrix.}
\item{\code{dist.eigenval }}{returns the eigenvalues of the dissimilarity
matrix.}
\item{\code{dissim.ord }}{returns a list: the first item is the the 
ordination distances, the second is the dissimilarity matrix distances.
}
\item{\code{dissim.GOF }}{returns the goodness of fit values of the
dissimilarity matrix, for various numbers of dimensions used.
}
\item{\code{dissim.tree }}{returns a list: the first item is the tree 
distances, the second is the dissimilarity matrix distances.
}
\item{\code{dissim.pvar }}{returns a numeric vector containing the 
percent variation explained by each axis (where each sample corresponds 
to an axis).
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson
}
%\note{
%}

\seealso{
\code{\link{vegdist}}, \code{\link{hclust}}, \code{\link{dissim.plot}}
}
\examples{
data(ITS1)

# calculate clustering, using default method
dissim.clust(ITS1)

# calculate tree distances, specifying a distance method
# (but use default clustering method)
dissim.tree(ITS1, dist.method="euclidean")

# calcualte ordination distances, specifying both distance 
# and ordination methods
dissim.ord(ITS1, dist.method="bray", k=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
