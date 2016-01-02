\name{combine.OTU}
\alias{combine.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine Non Overlapped OTU tables From The Same Community
}
\description{
This function combines otu tables from the same community but based
on independent sequencing runs. Such combined otu table gives a
more complete profile of the microbial community than each individual
otu table does.
This function should NOT be used to combine ITS1 and ITS2 otu
tables if they were extracted from long NGS sequences.
}
\usage{
combine.OTU(data, meta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list of otu tables to be combined.
}
  \item{meta}{
the metadata that should have the same number and order of the samples
as the otu tables do.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{combine.OTU} returns a data frame of combined otu tables which
have the same samples. Samples in the output will match those in
the metadata provided.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
%\note{
%
%}
%% ~Make other sections like Warning with \section{Warning }
%% {....} ~
\seealso{
\code{\link{match.data}}
}
\examples{
data(ITS1, ITS2, meta)
meta.new <- head(meta)
\dontrun{
# for demonstration purposes only, Not recommend to combine
# ITS1 and ITS2 otu tables that both regions were extracted from
# long NGS sequences
comb <- combine.OTU(data=list(ITS1=ITS1, ITS2=ITS2), meta=meta.new)
stopifnot(identical(colnames(comb)[1:(ncol(comb)-1)],
                    rownames(meta.new)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
