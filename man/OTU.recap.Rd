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
OTU.recap(data, ranks=c("p", "c", "o", "f", "g"), 
          brewer.pal="Pastel1", file=NULL, ext="pdf", 
          width=12, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of otu tables. See also \code{\link{RAM.input.formatting}}.
}
  \item{ranks}{
a vector of taxonomic ranks. See also \code{\link{RAM.rank.formatting}}
}
  \item{brewer.pal}{
one of the color patterns available in RColorBrewer. 
See \code{\link{brewer.pal}} for available selections.
}
  \item{file}{
filename to save the plot.
}
  \item{ext}{
extention of the filename to save the plot.
}
  \item{width}{
width of the plot
}
  \item{height}{
heigth of the plot
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function returns either a data frame or a list of data 
frames. If a single otu was provided, it returns the a dataframe
with information of how many otuIDs and sequences being classified 
at selected taxonomic ranks. If more than 1 otu tables being provided,
it returns a list, with the first a few are data frames of 
classification summary of each otu table, the last is a list showing 
taxa found only in one of the otu data set.

This function also generates a barplot for the percent classified
otus and sequences at each given rank.
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
\code{\link{RAM.rank.formatting}} and \code{\link{RAM.input.formatting}}.
}
\examples{
data(ITS1, ITS2)

ranks <- c("p", "c", "o", "f", "g")
df <- OTU.recap(data=list(ITS1=ITS1, ITS2=ITS2), ranks=ranks)
class(df)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
