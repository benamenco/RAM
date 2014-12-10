\name{filter.META}
\alias{filter.META}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select METADATA Variables
}
\description{
This function will remove metadata variables with only one level 
and /or remove variables with missing data or neither numeric nor 
factor/character (NNF). 
}
\usage{
filter.META(meta=meta, excl.na=TRUE, excl.NNF=TRUE, exclude=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta}{
the metadata table to be analyzed.
}
  \item{excl.na}{
logical. Whether or not remove variables that contain missing data. 
}
  \item{excl.NNF}{
logical. Whether or not remove variables that neither are numeric 
nor factor/character. 
}
  \item{exclude}{
optional. If is \code{NULL}, the function only removes variables 
with only one level or NNF. Otherwise, the variables in the 
\code{exclude} will also be removed from the metadata table.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a data frame with the 
following metadata variables being removed: 1) with missing data; 
2) NNF if \code{excl.NNF} is \code{TRUE}; and 3) in the 
\code{exclude} list. 
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
data(meta)
\dontrun{  
# add a new column with NA 
meta.nw <- meta
meta.nw$na <- c(rep(NA, nrow(meta.nw)-3), c(1, 3, 5))
meta.nw$nf <- rep("Canada", nrow(meta.nw))

meta.fil <- filter.META(meta.nw)
meta.fil <- filter.META(meta.nw, excl.na=FALSE, excl.NNF=FALSE, 
                        exclude=c("Province", "Latitude"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
