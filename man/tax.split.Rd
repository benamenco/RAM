\name{tax.split}
\alias{tax.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split OTU Tables By Taxonomic Rank
}
\description{
This function consumes an OTU table and splits its taxonomy 
columns into the seven major taxonomic ranks. It returns a data 
frame preserving all numerical data, but changing the 'taxonomy'
column to the name of the appropriate rank, and preserving only 
the classifications at that rank.
}
\usage{
tax.split(otu1, otu2 = NULL, rank = NULL)
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
the (optional) rank to split at and return 
(see ?RAM.rank.formatting for formatting details).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
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
        frames (otu1 and otu2 split at the given rank).
    \item If \code{rank} is not given: a list containing two 
        lists. The first sublist represents \code{otu1}, the 
        second \code{otu2}. The sublists contain seven data 
        frames, which are the data split at each taxonomic rank 
        (see Examples).
  }
\item If \code{otu2} is not given:
  \itemize{
    \item If \code{rank} is given: a single data frame 
       (otu1 split at the given rank).
    \item If \code{rank} is not given: a list containing seven 
       data frames (otu1 split at each taxonomic rank).
  }
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
\note{
This function may seem similar to \code{get.rank}, but they are
distinct. \code{get.rank} only returns the entries classified 
at that taxonomic rank, and so some OTUs might be omitted in 
the returned data frame. With \code{tax.split}, it is guaranteed
that all OTUs will be preserved in the returned data table 
(although they may be missing taxonomic classification at that 
rank). If no OTUs are classified at the given rank, the taxonomy
column for that rank will be filled with empty strings.
}
\seealso{
\code{\link{get.rank}}
}
\examples{
data(ITS1, ITS2)
# split only ITS1 data at a single rank
ITS1.tax.p <- tax.split(ITS1, rank="phylum")
# split only ITS1 data at all ranks
ITS1.tax.all <- tax.split(ITS1)
# split ITS1 and IST2 data at a given rank
lst <- tax.split(ITS1, ITS2, rank="c")
# split ITS1 and ITS2 at every rank
lst.all <- tax.split(ITS1, ITS2)
stopifnot(identical(lst.all$otu1$phylum, ITS1.tax.p))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
