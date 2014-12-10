\name{fread.OTU}
\alias{fread.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Load Large OTU Table
}
\description{
This function is same as \code{\link{read.OTU}} except using 
\code{\link{fread}} for loading large data sets.
}
\usage{
fread.OTU(file, sep="auto")
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
Wen Chen.
}
\note{
The OTU table should only contain classifications for the seven 
major taxonomic ranks, additional ranks will break some functions 
in the package. To remove those other classifications, try running 

\code{sed -i.backup -e 's/s[a-z]__[^;]*; //g' -e 's/t__[^;]*; //g' $FILE}

where \code{$FILE} is your OTU table. The letter \code{t} can be 
replaced in the second expression for any other letter which appears 
as a prefix. For example, adding \code{-e 's/u__[^;]*; //g'} 
before \code{$FILE} would remove any entries formatted like 
\code{u__test_classification; }.

Additionally, if your OTU table starts with the entry \code{#OTU ID}, 
that cell needs to be removed before the table can be read with 
\code{fread.OTU}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.OTU}}, \code{\link{fread}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
    my.OTU <- fread.OTU("path/to/otu")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
\keyword{IO}
