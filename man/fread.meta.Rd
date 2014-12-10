\name{fread.meta}
\alias{fread.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load Metadata Table
}
\description{
This function is same as \code{\link{read.meta}} to read in 
data; except using \code{\link{fread}} for loading large data sets.
}
\usage{
fread.meta(file, sep="auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
a character vector specifying the file path to your file.
}
  \item{sep}{
the character used to separate cells in the file.  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a data frame with the information from the file. The first 
row and column are used for the names of the other rows and columns.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.meta}}, \code{\link{fread}}
}
\examples{
  \dontrun{
    my.meta <- fread.meta("path/to/meta")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
\keyword{IO}
