\name{correlation}
\alias{correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Of Correlation Coefficient 
}
\description{
This function plot correlation relationship among taxa at a give 
rank and / or numeric variables of metadata.
}
\usage{
correlation(data=NULL, is.OTU=TRUE, meta=NULL, rank="g", 
            sel=NULL, sel.OTU=TRUE, data.trans=NULL, 
            method="pearson", main=NULL, file=NULL, 
            width=8, height=8)   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame that either an OTU table or taxonomy abundance 
matrix, can be missing but if metadata is also missing, an error 
message will
be raised.
}
  \item{is.OTU}{
logical. Whether or not the data is an OTU table.
}
  \item{meta}{
the metadata table to be used.
}
  \item{rank}{
the taxonomic rank to use (see ?RAM.rank.formatting for 
formatting details). 
}
  \item{sel}{
optional. It is a character vector of selected otuIDs or taxa 
names at a given taxonomic rank. If provided, \code{sel.OTU} 
should be set to decribe the type of IDs, i.e. \code{TRUE} means 
otuIDs, \code{FALSE} means taxa names. If provide, only the 
selected taxa will be ploted; otherwise, all taxa will be ploted.
}
  \item{sel.OTU}{
logical. Whether or not the selected items from data are otuIDs. 
If \code{FALSE}, \code{sel} should be a string vector of taxa 
names at a given \code{rank}.
}
  \item{data.trans}{
a character string of one of the following, "total", "log", 
"hellinger" etc, see \code{?vegan::decostand} for details and 
other data transformation methods.
}
  \item{method}{
a character string, can be one of the following, "pearson", 
"kendall", "spearman" for the calculation of correlation 
coefficient (or covariance) is to be computed (see 
\code{?stats::cor} for details)  
}
  \item{main}{
a character string. The title of the plot.  
}
  \item{file}{
the file path where the image should be created 
(see ?RAM.plotting).  
}
  \item{height}{
the height of the image to be created (in inches).  
}
  \item{width}{
the width of the image to be created (in inches).   
}
}
\details{
This function uses \code{stats::cor} to calculate correlation 
coefficient (or covariance), and uses \code{lattice::levelplot} 
to generate the graph. (see References)

Option \code{sel} is optional, however, it raises an error if 
the total number of variables to be plotted was too big, and no 
plot will be generated. 
}
\value{
This function generates a graph showing correlation relationship 
among OTUs or taxa at a given rank, and numeric variables of 
metadata
}
\references{
Sarkar, Deepayan (2008) _Lattice: Multivariate Data Visualization
with R_, Springer.  <URL: http://lmdvr.r-forge.r-project.org/>

Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) _The New S
Language_.  Wadsworth & Brooks/Cole.
  
}
\author{
Wen Chen.
}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link[stats]{cor}}
\code{\link[lattice]{levelplot}}
}
\examples{
data(ITS1, meta)
    
# only plot the first 10 OTUs
sel <- rownames(ITS1)[1:10]
correlation(data=ITS1, meta=meta, is.OTU=TRUE, sel.OTU=TRUE, 
            sel=sel, rank="g")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
