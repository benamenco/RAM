\name{assist.ado}
\alias{assist.ado}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform ADONIS Analysis for OTU Tables Or Taxonomic Abundance 
Matrix
}
\description{
This function simplifies ADONIS analysis by abstracting away 
some of the complexity and returning a list of useful measures.
}
\usage{
assist.ado(data, meta, is.OTU=TRUE, ranks=NULL,           
           data.trans=NULL, dist=NULL, meta.strata=NULL, 
           perm=1000, top=NULL, mode="number")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an OTU table or a taxonomy abundance matrix.
}
  \item{is.OTU}{
logical. If the data is an OTU table, set \code{is.OTU} \code{TRUE};
otherwise, set it as \code{FALSE}.
}
  \item{meta}{
the metadata table to be used (must have same samples as 
\code{data}.
}
  \item{ranks}{
optional. If \code{ranks} is not provided, will test for OTUs, 
otherwise, will test on taxa at defined ranks. If \code{data} is 
a taxonomic abundance matrix, \code{ranks} can be \code{NULL}
}
  \item{data.trans}{
optional. Transform the \code{data} using method from the function 
\code{\link{decostand}} 
}
  \item{dist}{
optional. the name of any method used in \code{\link{vegdist}} to 
calculate pairwise distances. See also \code{\link{adonis}} and 
\code{\link{vegdist}}.
}
  \item{meta.strata}{
optional. A metadata variable within which to constrain permutations.
See also \code{\link{adonis}}
}
  \item{perm}{
a numeric number of replicate permutations used for the hypothesis 
test used in \code{\link{adonis}}. 
}
  \item{top}{
optional. Select the top taxa or OTUs. 
See also \code{\link{data.revamp}}
}
  \item{mode}{
a character vector, one of "\code{percent}" or "\code{number}". If
\code{number}, then \code{top} groups will be selected based on total
sequence count. If \code{percent}, then top groups will be selected
based on relative abundance. See also \code{\link{data.revamp}}
}
}
%\details{
%%  Maybe explain full/restricted models
%}
\value{
This function returns a list containing outputs from 
\code{\link{adonis}} test. 
\itemize{
\item If \code{is.OTU} is \code{TRUE} and \code{ranks} is not given: 
      the output is a length one list named \code{LCA_OTU}. 
\item If \code{is.OTU} is \code{TRUE} and \code{ranks} is given:
      the output is a list with a length same as the number of 
      taxonomic ranks provided. Each member of the list is named 
      after the rank it processed at. 
\item If \code{is.OTU} is \code{FALSE}, the output is a length 
      one list named \code{Taxa}. 
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{adonis}}
}
\examples{
data(ITS1, meta)
\dontrun{
# test OTUs
data <- list(ITS1=ITS1, ITS2=ITS2)
assist.ado(data, is.OTU=TRUE,meta=meta, ranks=NULL, 
           data.trans="log", dist=NULL)

# test taxa at different ranks
ranks <- c("p", "c", "o", "f", "g")
ado <- assist.ado(data, is.OTU=TRUE, meta=meta, ranks=ranks, 
           data.trans="log", dist="bray" )

# test genera
g1 <- tax.abund(data, rank="g", drop.unclassified=TRUE)
data <- list(g1=g1)
assist.ado(data, is.OTU=FALSE, meta=meta, ranks=NULL, 
           data.trans="log", dist="bray" )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
