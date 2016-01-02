\name{RAM-package}
\alias{RAM-package}
\alias{RAM}
\docType{package}
\title{
Analysis of Amplicon-Based Metagenomic Data
}
\description{
The RAM package provides a series of functions to make amplicon
based metagenomic analysis more accessible. The package is 
designed especially for those who have little or no experience 
with R. This package calls heavily upon other packages (such as 
\code{vegan} and \code{ggplot2}), but the functions in this 
package either extend their functionality, or increase the 
ease-of-use.
}
\details{
\tabular{ll}{
Package: \tab RAM\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2014-12-10\cr
License: \tab MIT License, Copyright (c) 2014 Government of 
Canada\cr
}
Load data from .csv-formatted OTU files with \code{read.OTU} or
\code{fread.OTU}, then process the data with other commands.
Type the command \code{library(help = RAM)} for a full index of
all help topics, or \code{ls("package:RAM")} to get a list of 
all functions in the package.
Type \code{data(ITS1, ITS2, meta)} to load sample data sets of 
RAM, which include the following data of 16 samples:
1) ITS1: OTU table of fungal internal transcribed spacer region1
3) ITS2: OTU table of fungal internal transcribed spacer region2
3) meta: associated metadata
4) alignment for \code{\link{seq_var}}
Type \code{citation("RAM")} for how to cite this package.
This pacakge contains information licensed under the Open 
Government Licence - Canada. See \code{\link{group.spatial}} 
for further details.
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
\code{\link[vegan]{vegdist}}, \code{\link[ggplot2]{ggplot}}
}
\examples{
\dontrun{
# load data from your own files...
otu1 <- fread.OTU("path/to/OTU/table")
otu2 <- read.OTU("path/to/OTU/table")
meta1 <- fread.meta("path/to/meta/table")
meta2 <- read.meta("path/to/meta/table")
# ...or use the included sample data
data(ITS1, ITS2, meta)
data <- list(ITS1=ITS1, ITS2=ITS2)
dissim.heatmap(ITS1, meta, row.factor=c(City="City"))
dissim.alleig.plot(data)
data(alignment)
# type library(help = RAM) to get a full listing of help 
documents
}
}
