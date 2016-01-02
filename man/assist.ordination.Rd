\name{assist.ordination}
\alias{assist.ordination}
\alias{assist.cca}
\alias{assist.rda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform CCA and RDA Analysis for OTU Tables
}
\description{
This function simplifies CCA and RDA analysis by abstracting 
away some of the complexity and returning a list of useful 
measures.
}
\usage{
assist.cca(otu1, otu2 = NULL, meta, full = TRUE, exclude = NULL,
           rank, na.action=na.exclude)
assist.rda(otu1, otu2 = NULL, meta, full = TRUE, exclude = NULL,
           rank, na.action=na.exclude)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be used.
}
  \item{otu2}{
the second OTU table to be used.
}
  \item{meta}{
the metadata table to be used (must have same samples as
\code{otu1}/\code{otu2}).
}
  \item{full}{
logical. Should a full model be considered?
(If not, a restricted model is used).
}
  \item{exclude}{
A vector, either numeric or logical, specifying the columns to 
be removed from meta. If a character vector, columns with those 
names will be removed; if a numeric vector, columns with those 
indices will be removed.
}
  \item{rank}{
a character vector representing a rank. Must be in one of three
specific formats (see ?RAM.rank.formatting for help).
}
  \item{na.action}{
choice of one of the following: "na.fail", "na.omit" or 
"na.exclude", see \code{na.action} in \code{\link[vegan]{cca}} 
for detail.
}
}
%\details{
%%  Maybe explain full/restricted models
%}
\value{
If both \code{otu1} and \code{otu2} are given, a list of length 
2 will be returned with the following items (if only \code{otu1}
is given, a list of length 1 will be returned with these items):
\item{$GOF }{the goodness of fit scores for the model.}
\item{$VIF }{the VIF scores for the model.}
\item{$percent_variation }{
  the percent variation explained by each axis
}
\item{$CCA_eig }{Eigenvalues for CCA axes.}
\item{$CA_eig }{Eigenvalues for CA axes.}
\item{$anova }{the ANOVA results for the model.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%%  ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link[vegan]{cca}}, \code{\link[vegan]{anova.cca}}
}
\examples{
data(ITS1, meta)
cca.help <- assist.cca(ITS1, meta=meta, rank="p")
cca.help$anova
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
