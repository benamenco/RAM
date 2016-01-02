\name{get.rank}
\alias{get.rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get OTUs Classified at Taxonomic Rank(s)
}
\description{
This function returns the OTUs of the given OTU table(s) which 
areclassified at the given taxonomic rank.
}
\usage{
get.rank(otu1, otu2 = NULL, rank = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be used.
}
  \item{otu2}{
the second OTU table to be used.
}
  \item{rank}{
a character vector representing a rank. Must be in one of three
specific formats (see \code{\link{RAM.rank.formatting}} for 
help). If omitted, the function will repeat for all seven major 
taxonomic ranks.
}
}
%\details{
%}
\value{
The value returned by this function may become nested lists, so 
for convenience, any nested lists have been given descriptive 
items names to make accessing its elements simple 
(see Examples).
\itemize{
\item If \code{otu2} is given:
  \itemize{
    \item If \code{rank} is given: a list containing two data 
          frames (otu1 and otu2 selected at the given rank).
    \item If \code{rank} is not given: a list containing two 
          lists. The first sublist represents \code{otu1}, the 
          second \code{otu2}. The sublists contain seven data 
          frames, which are the OTU tables selected at each 
          taxonomic rank (see Examples).
  }
\item If \code{otu2} is not given:
  \itemize{
    \item If \code{rank} is given: a single data frame (otu1 
          selected at the given rank).
    \item If \code{rank} is not given: a list containing seven 
          data frames (otu1 selected at each taxonomic rank).
  }
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%%  ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
%\seealso{
%
%}
\examples{
    data(ITS1, ITS2)
    # the following are equivalent:
    ITS1.p <- get.rank(ITS1, rank="p")
    # this list has get.rank(ITS1, rank="k"),
    #               get.rank(ITS1, rank="p"), ...
    lst <- get.rank(ITS1)
    stopifnot(identical(ITS1.p, lst$phylum))
    # true
    # get a list of length 2: the item holds all ITS1 data, the
    # second holds ITS2 data
    lst.all <- get.rank(ITS1, ITS2)
    stopifnot(identical(ITS1.p, lst.all$otu1$phylum))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
