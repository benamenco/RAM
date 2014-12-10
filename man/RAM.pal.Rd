\name{RAM.pal}
\alias{RAM.pal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creat Color Palette
}
\description{
This function creates color palette, especially if the number of 
colors required is more than 12.
}
\usage{
RAM.pal(cols.needed=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cols.needed}{
an integer.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
col <- RAM.pal(40)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
