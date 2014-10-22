\name{OTU.recap}
\alias{OTU.recap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize OTU
}
\description{
This function summarize OTU table at each given taxonomic ranks. 
}
\usage{
OTU.recap(otu1, otu2=NULL, ranks="", labels=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be analyzed.
}
  \item{otu2}{
optional. the second OTU table to be used.
}
  \item{ranks}{
a vector of taxonomic ranks.
}
  \item{labels}{
gives labels to each OTU table.  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function returns either a data frame or a list of data 
frames
\itemize{
\item If \code{otu2} is given: a list containing the following 
          members: 1) summary of otu1; 2) summary of otu2, 
          3) taxa only in one OTU table.
\item If \code{otu2} is not given: a data frame summarizes otu1.
  }
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
\note{
 warning is raised when run strsplit() and can be ignored.
}

%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link{RAM.rank.formatting}}
}
\examples{
data(ITS1, ITS2)

df <- OTU.recap(otu1=ITS1, otu2=ITS2, 
                ranks=c("p","c","o","f", "g"), 
                labels=c("ITS1", "ITS2"))
class(df)
length(df)

df<-OTU.recap(ITS1, ranks=c("p","c","o","f", "g"), 
              labels=c("ITS1"))
class(df)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
