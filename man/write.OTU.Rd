\name{write.OTU}
\alias{write.OTU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write OTU-Style Table to .csv File 
}
\description{
Creates a .csv-formatted file with the data from the specified OTU 
table \code{data}. The file will be created in your current working 
directory (see ?getwd and ?setwd), unless specified otherwise by 
\code{file}. Note that if the \code{file} field does not end in 
.csv, ".csv" will be appended to the end of \code{file}.
}
\usage{
write.OTU(data, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The OTU table containing the data to be used.
}
  \item{file}{
The name of the .csv file to be created.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
  
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{write.csv}}, \code{\link{getwd}}, \code{\link{setwd}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(ITS1)
\dontrun{
write.OTU(ITS1, "my_file_name.csv")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}% __ONLY ONE__ keyword per line
