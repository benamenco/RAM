\name{shared.OTU}
\alias{shared.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of Shared OTUs Across ALL Subjects
}
\description{
This function consumes OTU tables and returns a list summarizing 
information about the presence of the OTUs in samples.
}
\usage{
shared.OTU(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of OTU tables to be analyzed.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{shared.OTU} returns a list containing the information 
calculated. The names associated with the list describe what that 
number represents; i.e. "#_of_OTUs_in_all_samples" shows how many 
OTUs in the given table were found to be present in all samples. 
The last item in the list is a character vector, containing the 
OTU number and taxonomic information of each OTU which was present 
in all samples. All entries in that column are of the form 
"OTU-taxonomic_classification".
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
\note{
The OTUs are determined to be absent/present using the "pa" method 
from the function \code{\link{decostand}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{decostand}}
}
\examples{
data(ITS1)
\dontrun{
shared <- shared.OTU(data=list(ITS1=ITS1))
shared <- shared.OTU(data=list(ITS1=ITS1, ITS2=ITS2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
