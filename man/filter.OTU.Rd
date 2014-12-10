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
filter.OTU(data, percent=NULL, number=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of otu tables to be processed.
See also \code{\link{RAM.input.formatting}}
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
The value returned by this function is a list of filtered otu 
tables provided by the user
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
data<-list(ITS1=ITS1, ITS2=ITS2)
\dontrun{
otu.001 <- filter.OTU(data=data, percent=0.01)
length(otu.001)
names(otu.001)
otu.50 <- filter.OTU(data=data, number=50)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
