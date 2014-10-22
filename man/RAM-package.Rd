\name{RAM-package}
\alias{RAM-package}
\alias{RAM}
\docType{package}
\title{
Analysis of Amplicon-Based Metagenomic Data
}
\description{
The RAM package provides a series of functions to make amplicon 
based metagenomic analysis more accessible. The package is designed 
especially for those who have little or no experience with R. This 
package calls heavily upon other packages (such as \code{vegan} and 
\code{ggplot2}), but the functions in this package either extend 
their functionality, or increase the ease-of-use. 
}
\details{
\tabular{ll}{
Package: \tab RAM\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2014-10-20\cr
License: \tab MIT License, Copyright (c) 2014 Government of Canada\cr
}
Load data from .csv-formatted OTU files with \code{read.OTU}, then 
process the data with other commands. Type the command 
\code{library(help = RAM)} for a full index of all help topics, or 
\code{ls("package:RAM")} to get a list of all functions in the package.

Type \code{data(ITS1, ITS2, meta)} to load a sample (fungal) data set.

Type \code{citation("RAM")} for how to cite this package.

This pacakge contains information licensed under the Open Government 
Licence - Canada. See \code{\link{group.spatial}} for further details.
}
\author{
Wen Chen and Joshua Simpson.

Maintainer: Wen Chen <wen.chen@agr.gc.ca>
}
%\references{
%~~ Literature or other references for background information ~~
%}

\keyword{package}
\seealso{
\code{\link[vegan:vegan-package]{vegan}},
\code{\link[ggplot2:ggplot2-package]{ggplot2}}
}
\examples{
\dontrun{
# load data from your own files...
ITS1 <- read.OTU("path/to/OTU/table")
ITS2 <- read.OTU("path/to/OTU/table")
meta <- read.meta("path/to/meta/table")}

# ...or use the included sample data 
data(ITS1, ITS2, meta)

dissim.heatmap(ITS1, meta, row.factor=c(City="City"))
dissim.alleig.plot(ITS1, ITS2)

# type library(help = RAM) to get a full listing of help documents
}
