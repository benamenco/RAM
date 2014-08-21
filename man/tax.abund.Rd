\name{tax.abund}
\alias{tax.abund}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregate OTU Data Based on Taxonomy
}
\description{
This function consumes OTU table(s) and (optionally)  a taxonomic rank, then
extracts the classification of each OTU at the given taxonomic rank, groups by
classification at the given rank, removes all groups with only 0 counts,
optionally removes all unclassified groups, sorts groups based on abundance, and
then returns the
transpose.
}
\usage{
tax.abund(otu1, otu2=NULL, rank=NULL, drop.unclassified=FALSE, 
          top=NULL, count=TRUE, mode="number")
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
a character vector representing a rank. Must be in one of three specific formats
(see ?RAM.rank.formatting for help). If omitted, the function will repeat for
all seven major taxonomic ranks.
}
  \item{drop.unclassified}{
logical. Should OTUs labelled "unclassified" or missing classification at the
given taxonomic rank be excluded?
}
  \item{top}{
the number of groups to select, starting with the most abundant. If \code{NULL},
all are selected.
}
  \item{count}{
logical. Should the actual count of each OTU be shown, or should the relative
abundances be shown?
}
  \item{mode}{
a character vector, one of "\code{percent}" or "\code{number}". If
\code{number}, then \code{top} many groups will be selected. If \code{percent},
then all groups with relative abundance in at least one sample above \code{top}
will be selected.

}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function may become nested lists, so for convenience,
any nested lists have been given descriptive items names to make accessing its
elements simple (see Examples).
\itemize{
\item If \code{otu2} is given:
  \itemize{
    \item If \code{rank} is given: a list containing two data frames (otu1 and
	  otu2 aggregated at the given rank).
    
    \item If \code{rank} is not given: a list containing two lists. The first
	  sublist represents \code{otu1}, the second \code{otu2}. The sublists
	  contain seven data frames, the aggregation of the data at each
	  taxonomic rank (see Examples).
		      }

\item If \code{otu2} is not given:
  \itemize{
    \item If \code{rank} is given: a single data frame (otu1 aggregated at the
	  given rank).
    
    \item If \code{rank} is not given: a list containing seven data frames (otu1
	  aggregated at each taxonomic rank).
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{RAM.rank.formatting}}
}
\examples{
data(ITS1, ITS2)
# aggregate based on phylum
ITS1.p <- tax.abund(ITS1, rank="p")

# aggregate based on all ranks; ignoring all unclassified OTUs
ITS1.taxa <- tax.abund(ITS1, drop.unclassified=FALSE)

# aggregate for one rank for both ITS1 and ITS2
lst <- tax.abund(ITS1, ITS2, rank="class")

# aggregate for all ranks for both ITS1 and ITS2
lst.all <- tax.abund(ITS1, ITS2)

stopifnot(identical(lst.all$otu1$phylum, ITS1.p))

# get the counts for all genera with relative abundance > 25\%
tax.abund(ITS1, rank="g", top=25, mode="percent", count=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
