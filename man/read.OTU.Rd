\name{read.OTU}
\alias{read.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Open OTU Table
}
\description{
Opens the given file and returns a data frame representing the
OTU table. This function use \code{\link{read.table}} function 
so is quite slow for large data sets, for which we recommend to 
use \code{\link{fread.OTU}} instead.
}
\usage{
read.OTU(file, sep=",")
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
Wen Chen and Joshua Simpson.
}
\note{
The OTU table should only contain classifications for the seven
major taxonomic ranks, additional ranks will break some 
functions in the package. To remove those other classifications,
try running \code{sed -i.backup -e 's/s[a-z]__[^;]*; 
//g' -e 's/t__[^;]*; //g' $FILE} where \code{$FILE} is your OTU 
table. The letter \code{t} can be replaced in the second 
expression for any other letter which appears as a prefix. For 
example, adding \code{-e 's/u__[^;]*; //g'} before \code{$FILE} 
would remove any entries formatted like 
\code{u__test_classification;}. Additionally, if your OTU table 
starts with the entry \code{#OTU ID}, that cell needs to be 
removed before the table can be read with \code{read.OTU}.
}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link{getwd}}, \code{\link{setwd}}, 
\code{\link{read.table}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
my.OTU <- read.OTU("path/to/otu", sep=",")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
\keyword{IO}
