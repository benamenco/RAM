\name{data.subset}
\alias{data.subset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset OTU And Metadata
}
\description{
This function subset OTUs and metadata based on user defined
values of metadata variables.
}
\usage{
data.subset(data, meta, factors="", values="", and=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of otu tables to be processed.
See also \code{\link{RAM.input.formatting}}.
}
  \item{meta}{
the metadata for subset.
}
  \item{factors}{
a vector containing metadata variables.
}
  \item{values}{
a vector containing values of interest in metadata variables.
}
  \item{and}{
logical. Determine whether all conditions needs to be met or not.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a list containing otu 
tables matching the filtering requirement. The last item in the 
output list is the associated new metadata table fit the 
requirement.
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
%\seealso{
% ~~further links~~
%}
\examples{
data(ITS1, ITS2, meta)
names(meta)
factors <- c("City", "Harvestmethod")
values <- c("City1", "Method1")
# match all requirements, and=TRUE
sub <- data.subset(data=list(ITS1=ITS1, ITS2=ITS2), meta=meta,
                   factors=factors, values=values, and=TRUE)
# match either of the requirements, and=FALSE
sub <- data.subset(data=list(ITS1=ITS1, ITS2=ITS2), meta=meta,
                   factors=factors, values=values, and=FALSE)
\dontrun{
names(sub)
ITS1.sub <- sub[["ITS1"]]
ITS2.sub <- sub[["ITS2"]]
meta.sub <- sub[["meta"]]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
