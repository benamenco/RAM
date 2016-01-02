\name{alignment}
\alias{alignment}
\docType{data}
\title{
Sample Alignment
}
\description{
This is an alignment for \code{\link{seq_var}}
package.
}
\usage{data(alignment)}
\format{
  An alignment with sequence ID being formatted as follows:
genus_name:accession:genus:species:strain_info/seqBegin-seqEnd.
The location of each party can be rearranged, and the separator 
can be other speciall characters, such as "|".
}
%\details{
%%~~If necessary, more details than the __description__ above~~
%}
\source{
Wen Chen
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(alignment)
str(alignment)
alignment
}
\keyword{datasets}
