\name{ITS1/ITS2}
\alias{ITS1/ITS2}
\alias{ITS1}
\alias{ITS2}
\docType{data}
\title{
Sample ITS1 and ITS2 Data
}
\description{
Sample ITS1 and ITS2 OTU tables.
}
\usage{
data(ITS1)
data(ITS2)
}
\format{
  A data frame with 4704 observations on the following 17 variables.
  \describe{
    \item{\code{P1001.1M1,  P1001.1M2,  P1001.1M3,  P1001.1M4,  P1001.1M5,  
                P1001.1M6,  P1001.1M7,  P1001.1M8,  P1001.1M9, P1001.1M10, 
                P1001.1M11, P1001.1M12, P1001.1M13, P1001.1M14, P1001.1M15, 
                P1001.1M16}}{Collection samples.}
    \item{\code{taxonomy}}{the taxonomic classification of the OTU.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Wen Chen, AAFC-AAC
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(ITS1, ITS2)

str(ITS1)
str(ITS2)
}
\keyword{datasets}
