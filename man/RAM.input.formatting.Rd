\name{RAM.input.formatting}
\alias{RAM.input.formatting}
\title{Data Input Formatting}
\description{
When use some \code{RAM} functions for the comparison of 
multiple OTU tables or taxonomic abundance matrices, the user 
needs to privde all input data sets as list with names being 
provided.
\describe{
  \item{one data set:}{data=list(data=otu)}
  \item{multiple data sets:}{
data=list(data1=otu1, data2=otu2, data3=otu3)}
  \item{an OTU table:}{
a data frame of otuIDs x sampleIDs with the last column named
"taxonomy"}
  \item{a taxonomy abundance matrix:}{
a data frame of sampleIDs x taxa (e.g. species)}
  \item{is.OTU:}{logical, many functions in RAM require the user
to set \code{is.OTU} to be \code{TRUE} for OTU tables or
\code{FALSE} for a taxonomy abundance matrices.}
}
}
\author{
Wen Chen.
}
%\seealso{
% --
%}
\examples{
data(ITS1, ITS2, meta)
# use otu tables
matched <- match.data(data=list(otu_ITS1=ITS1, otu_ITS2=ITS2),
           is.OTU=TRUE, meta=meta)
# taxonomy abundance matrices
g1 <- tax.abund(ITS1, rank="g")
g2 <- tax.abund(ITS2, rank="g")
matched <- match.data(data=list(genus_ITS1=g1, genus_ITS2=g2),
           is.OTU=FALSE, meta=meta)
}
