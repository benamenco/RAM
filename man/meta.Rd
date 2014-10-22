\name{meta}
\alias{meta}
\docType{data}
\title{
Sample Metadata for ITS1/ITS2
}
\description{
This data frame provides sample metadata for the ITS1/ITS2 data 
included in this
package.
}
\usage{data(meta)}
\format{
  A data frame with 16 observations on the following 10 variables.
  \describe{
    \item{\code{Sample}}{a factor with levels \code{Sample1} 
      \code{Sample2} \code{Sample3} \code{Sample4} \code{Sample5} 
      \code{Sample6} \code{Sample7} \code{Sample8}}
    \item{\code{City}}{a factor with levels \code{City1} \code{City2}}
    \item{\code{Crop}}{a factor with levels \code{Crop1} \code{Crop2}}
    \item{\code{Plots}}{a factor with levels \code{1} \code{2}}
    \item{\code{Harvestmethod}}{a factor with levels \code{Method1}
      \code{Method2}}
    \item{\code{Harvestdate}}{a Date}
    \item{\code{Ergosterol_ppm}}{a numeric vector}
    \item{\code{Province}}{a character vector}
    \item{\code{Latitude}}{a numeric vector}
    \item{\code{Longitude}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Wen Chen and Joshua Simpson.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(meta)

str(meta)
}
\keyword{datasets}
