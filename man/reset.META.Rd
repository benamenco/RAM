\name{reset.META}
\alias{reset.META}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reset  OTU
}
\description{
This function reset data type of metadata variables. 
}
\usage{
reset.META(meta, factor=NULL, numeric=NULL, date=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta}{
data frame. The metadata table to reset variable data type.
}
  \item{factor}{
a string or charactor vector, containing the column names of 
metadata variables to be set as factor.
}reset.META
  \item{numeric}{
a string or charactor vector, containing the column names of 
metadata variables to be set as numeric. 
}
  \item{date}{
a string or charactor vector containing the column names of 
metadata variables to be set as date.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function returns the same metadata with variables being reset 
to desired data type. Warnings or errors may be raise if the format 
of the original data cannot be recognized by R.  
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~further links~~
%}
\examples{
    data(meta)
    str(meta)

    # for demonstration purpose only
    meta.new <- reset.META(meta, factor=c("Plots"), 
                           numeric=c("City", "Province"))
    str(meta.new)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
