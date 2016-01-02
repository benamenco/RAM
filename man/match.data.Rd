\name{match.data}
\alias{match.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Match Samples In Ecology Data Sets and Metadata
}
\description{
This function will match samples in ecology data sets, either 
OTU tables or taxonomy abundance matrices, and those in 
metadata. It makes sure that datasets contains same samples in 
the same order.
}
\usage{
match.data(data, is.OTU=TRUE, meta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of ecology data sets. if \code{is.OTU} is \code{TRUE},
they should be OTU tables, otherwise should be taxonomy 
abundance matrices. See also \code{\link{RAM.input.formatting}}.
}
  \item{is.OTU}{
logical, whether or not the ecology data sets are OTU tables.
}
  \item{meta}{
metadata associated with input ecology data sets.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
%\note{
%%  ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link{RAM.input.formatting}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
data(ITS1, ITS2, meta)
meta <- meta[1:8, ]
# use otu tables
matched <- match.data(data=list(otu_ITS1=ITS1, otu_ITS2=ITS2),
           is.OTU=TRUE, meta=meta)
# use taxonomy abundance matrices
g1 <- tax.abund(ITS1, rank="g")
g2 <- tax.abund(ITS2, rank="g")
matched <- match.data(data=list(genus_ITS1=g1, genus_ITS2=g2),
           is.OTU=FALSE, meta=meta)
# class(matched)
# names(matched)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}% __ONLY ONE__ keyword per line
