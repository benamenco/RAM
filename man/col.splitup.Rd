\name{col.splitup}
\alias{col.splitup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split Column Of Data Frame
}
\description{
This function output consumes a data frame and split one by 
defined separator.
}
\usage{
col.splitup(df, col="", sep="", max=NULL, names=NULL, drop=TRUE)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
a data frame.  
}
  \item{col}{
name of a column in df.  
}
  \item{sep}{
the separator to split the column. It can be regular expression. 
}
  \item{max}{
optional. The number of columns to be split to.  
}
  \item{names}{
optional. The names for the new columns.  
}
  \item{drop}{
logical. Whether or not to keep the original column to be split 
in the output.  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a data frame. The selected 
column is split each separator and appended to the original data 
frame. The original column may may not to be kept in the output 
as defined by option \code{drop}.

The number of columns to be split to depends on three factors, 
1) the maximum columns that the original column can be split to 
by each separator; 2) the user definde \code{max}; and 3) the 
length of the column names defined by \code{names}. This function 
will split the column to the maximun number of the 3, empty 
columns will be filled with empty strings.  
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.  
}
%\note{
%
%}

%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

%\seealso{
%
%}
\examples{
data(ITS1)
    
# filter.OTU() returns a list
otu <- filter.OTU(ITS1, percent=0.001)[[1]]
# split and keep taxonomy column
otup.split <- col.splitup(otu, col="taxonomy", sep="; ", 
                          drop=FALSE)

# give new column names
tax.classes <- c("kingdom", "phylum", "class", 
                  "order", "family", "genus")
otup.split <- col.splitup(otu, col="taxonomy", sep="; ", 
                           drop=TRUE, names=tax.classes)
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
