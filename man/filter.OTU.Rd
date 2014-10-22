\name{filter.OTU}
\alias{filter.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter OTU
}
\description{
This function filter OTU table by counts or relative abundance. If 
filter by counts, otus having total counts more than a threshhold 
will be kept. If filter by relative abundance, otus with the maximum 
relative abundance greater than a threshhold in at least one subject 
will be kept. 
}
\usage{
filter.OTU(otu1, otu2=NULL, percent=NULL, number=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be analyzed.
}
  \item{otu2}{
optional. the second OTU table to be used.
}
  \item{percent}{
a floating point greater than 0 and less or equals to 1. 
Filter OTU table by relative abundance.   
}
  \item{number}{
an integer. Filter OTU table by counts.  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function may become nested lists, so 
for convenience, any nested lists have been given descriptive 
items names to make accessing its elements simple (see Examples).
\itemize{
\item If \code{otu2} is given: a list containing the following 
      members: 1) filtered otu1; 2) filtered otu2, 3) a list of 
      otuIDs of otu1 that kept; 4) a list of otuIDs of otu2 that 
      kept.
\item If \code{otu2} is not given: a list containing the following 
      members: 1) filtered otu1; 2) a list of otuIDs of otu1 that 
      kept.
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
    data(ITS1, ITS2, meta)

    otu.fil <- filter.OTU(ITS1, otu2=ITS2, percent=0.01)
    otu.fil <- filter.OTU(ITS1, otu2=ITS2, number=50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
