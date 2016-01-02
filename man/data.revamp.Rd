\name{data.revamp}
\alias{data.revamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform OTU Table
}
\description{
This function consumes and transforms either an OTU table or a
taxonomy abundance matrix. If an OTU table was provided, it will
be either transposed without the "taxonomy" column, but each 
otuID will be renames with it's LCA classification appended; 
or being transformed to be taxonomic abundance matrix at the 
ranks set by \code{ranks}.
If a taxonomic abundance matrix is provided, it will be kept the
same with proper data transformation as defined by 
\code{stand.method} option.
}
\usage{
data.revamp(data, is.OTU=TRUE, ranks=NULL, stand.method=NULL,
            top=NULL, mode="number")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an OTU table or a taxonomic abundance matrix.
}
  \item{is.OTU}{
logical. If an OTU table was provided, \code{is.OTU} should be 
set as \code{TRUE}; otherwise, it should be set as \code{FALSE}.
}
  \item{ranks}{
optional. If no \code{ranks} was provided, the OTU table will be
processed by \code{\link{LCA.OTU}} and then transposed with 
sampleIDs being row names and otuIDs being column names. If 
\code{ranks} was provided, the OTU table will be processed by 
\code{\link{tax.abund}} at each given taxonomic ranks.
See also \code{\link{RAM.rank.formatting}}. The unclassified 
taxon groups are removed.
}
  \item{stand.method}{
optional. Transform the output using method from the function
\code{\link{decostand}}
}
  \item{top}{
optional. Select the top taxa or OTUs.
}
  \item{mode}{
a character vector, one of "\code{percent}" or "\code{number}". 
If \code{number}, then \code{top} many groups will be selected. 
If \code{percent}, then all groups with relative abundance in 
at least one sample above \code{top} will be selected.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a list, so for 
convenience, any nested lists have been given descriptive items 
names to make accessing its elements simple (see Examples).
\itemize{
\item If \code{is.OTU} is \code{TRUE} and \code{ranks} is not 
      given: the output is a length one list named 
      \code{LCA_OTU}.
\item If \code{is.OTU} is \code{TRUE} and \code{ranks} is given:
      the output is a list with a length same as the number of
      taxonomic ranks provided. Each member of the list is named
      after the rank it processed at.
\item If \code{is.OTU} is \code{FALSE}, the output is a length
      one list named \code{Taxa}.
  }
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
% \note{
%% ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
%\seealso{
% ~~further links~~
%}
\examples{
data(ITS1, ITS2, meta)
data.new <- data.revamp(data=list(ITS1=ITS1), is.OTU=TRUE,
                        ranks=c("f", "g"), stand.method="log")
\dontrun{
data.new <- data.revamp(data=list(ITS1=ITS1), is.OTU=TRUE,
                        ranks=NULL, stand.method="log")
data.new <- data.revamp(data=list(ITS1=ITS1, ITS2=ITS2),
           is.OTU=TRUE, ranks=c("f", "g"), stand.method="total")
names(data.new)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
