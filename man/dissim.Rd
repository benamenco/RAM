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
These functions calculate different measures related to 
dissimilarity matrices. All of these functions allow you to 
specify one of many dissimilarity indices to be used.
}
\usage{
dissim.clust(elem, is.OTU=TRUE, stand.method=NULL,
             dist.method="morisita", clust.method="average")
dissim.eig(elem, is.OTU=TRUE, stand.method=NULL,
           dist.method="morisita")
dissim.ord(elem, is.OTU=TRUE, stand.method=NULL,
           dist.method="morisita", k=NULL)
dissim.GOF(elem, is.OTU=TRUE, stand.method=NULL,
           dist.method="morisita")
dissim.tree(elem, is.OTU=TRUE, stand.method=NULL,
           dist.method="morisita", clust.method="average")
dissim.pvar(elem, is.OTU=TRUE, stand.method=NULL,
            dist.method="morisita")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elem}{
an ecology data set that can be an OTU table or a taxonomy
abundance table. See \code{\link{RAM.input.formatting}} for 
details.
}
  \item{is.OTU}{
logical, whether the ecology data sets are OTU tables or 
taxonomy abundance matrices. 
See \code{\link{RAM.input.formatting}} for details.
}
  \item{stand.method}{
optional, if \code{is.null}, the standardization method for data
transforamtion; must be one of the following: "total", "max",
"frequency", "normalize", "range", "standardize", "pa",
"chi.square", "hellinger", "log".
See also \code{\link{decostand}}.
}
  \item{dist.method}{
the dissimilarity index to be used; one of "\code{manhattan}",
"\code{euclidean}", "\code{canberra}", "\code{bray}",
"\code{kulczynski}", "\code{jaccard}", "\code{gower}",
"\code{altGower}","\code{morisita}", "\code{horn}",
"\code{mountford}", "\code{raup}", "\code{binomial}",
"\code{chao}", or "\code{cao}".
See also \code{\link{vegdist}}.
}
  \item{k}{
the number of dimensions desired. If \code{NULL}, the maximum 
value will be calculated and used.
}
  \item{clust.method}{
the method used for clustering the data. Must be one of 
"\code{ward}", "\code{single}", "\code{complete}", 
"\code{average}", "\code{mcquitty}", "\code{median}", 
or "\code{centroid}". See also \code{\link{hclust}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{\code{dissim.clust }}{returns a hierarchical clustering 
of the dissimilarity matrix.}
\item{\code{dist.eigenval }}{returns the eigenvalues of the 
dissimilarity matrix.}
\item{\code{dissim.ord }}{returns a list: the first item is the 
the ordination distances, the second is the dissimilarity matrix
distances.
}
\item{\code{dissim.GOF }}{returns the goodness of fit values 
of the dissimilarity matrix, for various numbers of dimensions 
used.
}
\item{\code{dissim.tree }}{returns a list: the first item is 
the tree distances, the second is the dissimilarity matrix 
distances.
}
\item{\code{dissim.pvar }}{returns a numeric vector containing 
the percent variation explained by each axis (where each sample 
corresponds to an axis).
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
\code{\link{decostand}}, \code{\link{vegdist}},
\code{\link{hclust}}, \code{\link{dissim.plot}}
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
