\name{assist.NB}
\alias{assist.NB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Negative Binomial Test For OTUID or Taxon
}
\description{
This function does negative binomial test for a given otuID or taxon
}
\usage{
assist.NB(data, meta, is.OTU=TRUE, rank=NULL, meta.factors=NULL,
          anov.fac=NULL, taxon="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an ecology data set to be analyzed.
}
  \item{meta}{
the metadata table to be analyzed.
}
  \item{is.OTU}{
logical. If an OTU table was provided, \code{is.OTU} should be set 
as \code{TRUE}; otherwise, it should be set as \code{FALSE}.
}

  \item{rank}{
optional. If no \code{rank} was provided, the data will be used as
it is, if \code{rank} is provided, if data is an OTU table, it will
be converted to taxonomic abundance matrix at the given rank, no 
change will be made for a data that has already been a taxonomic 
abundance matrix. 
See also \code{\link{tax.abund}} and \code{\link{data.revamp}}
}
  \item{meta.factors}{
optional. If provided, will only test the model on selected metadata
variables; otherwise, will test all variables in the metadata table.
}
  \item{anov.fac}{
optional. Whether or not to do anova test on a metadata variable.
}
  \item{taxon}{
a length one charactor. Can either be an otuID or a taxon name.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function return a list of outputs of the negative bionomial
modeling for a selected otuID or taxa. Members of this output
list are: "NB.model", "tax.met", "taxon", "factors", "anova".
  \item{NB.model}{is the negative bionomial model}
  \item{tax.met}{is a dataframe with combined the taxon and metadata}
  \item{taxon}{is either a taxon name or in LCA_otuID format, 
see also \code{\link{LCA.OTU}}
}
  \item{factors}{shows which metadata variable had significant impact}
  \item{anova}{shows anova test of a metadata variable, this will not
be available if \code{anov.fac} is \code{NULL}
}
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
m <- meta[, c(2,3,5,7)]
\dontrun{  
# for usage demonstration purpose only, may not fit the negative 
# binomial distribution model.
nb <- assist.NB(ITS1, meta=m, rank="g", anov.fac="Harvestmethod", 
              taxon=rownames(ITS1)[1])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
