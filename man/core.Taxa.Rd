\name{core.Taxa}
\alias{core.Taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show Summary of Core Taxa
}
\description{
This function returns a list showing taxa at the given taxonomic
rank that present in a pre-defined percent of samples in each
level of a given metadata category.
}
\usage{
core.Taxa(data, is.OTU=FALSE, rank="g",
         drop.unclassified=TRUE,
         meta, meta.factor="", percent=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables or taxonomy abundance matrices.
}
  \item{is.OTU}{
logical. If \code{TRUE}, data is an OTU table;
otherwise a taxonomy abundance matrix should be provided.
}
  \item{rank}{
the taxonomic rank of classification (see ?RAM.rank.formatting 
for formatting details).
}
  \item{drop.unclassified}{
logical, whether or not exclude unclassified groups.
See also \code{\link{tax.abund}}
}
  \item{meta}{
the metadata table to be analyzed.
}
  \item{meta.factor}{
the metadata qualitative variable
}
  \item{percent}{
the percent of samples in each level of the given metadata 
variable
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{core.Taxa} returns a list containing taxa at a given rank
that present in a pre-defined percent of samples in each level 
of a given metadata category. The outputs describe the following
information for each level of a given metadata variable:
1) core taxa; 2) percent of core taxa sequences vs. total 
sequences in each levels of the given metadata variable.  The 
last item in the list show the same information of taxa that in 
all levels.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
\note{
The taxa are determined to be absent/present using the "pa" 
method from the function \code{decostand}.
}
\seealso{
\code{\link{decostand}}
}
\examples{
data(ITS1, meta)
# taxa shared by 50 percent samples of each city
core <- core.Taxa(data=list(ITS1=ITS1), is.OTU=TRUE, meta=meta,
                  rank="g", meta.factor="City", percent=0.5)
\dontrun{
data(ITS1, ITS2, meta)
core <- core.Taxa(data=list(ITS1=ITS1, ITS2=ITS2), is.OTU=TRUE,
                  meta=meta, rank="g", meta.factor="City",
                  percent=0.7)
# use taxonomy abundance matrix
g1<-tax.abund(ITS1, rank="g")
core <- core.Taxa(data=list(genus_ITS1=g1), is.OTU=FALSE, 
                  meta=meta, rank="g", meta.factor="City", 
                  percent=0.9)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
