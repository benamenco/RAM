\name{read.meta}
\alias{read.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Open Metadata Table
}
\description{
Opens the given file and return a data frame representing the 
metadata. This function use \code{\link{read.table}} to read in
data; for large data sets, we recommend \code{\link{read.meta}}.
}
\usage{
read.meta(file, sep=",")
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
Returns a data frame with the information from the file. The 
first row and column are used for the names of the other rows 
and columns.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson
}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link{read.meta}}, \code{\link{read.table}}
}
\examples{
  \dontrun{
    my.meta <- read.meta("path/to/meta")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
\keyword{IO}
