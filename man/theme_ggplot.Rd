\name{theme_ggplot}
\alias{RAM.border}
\alias{RAM.color}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Customized Themes For GGPLOT
}
\description{
RAM provides some customized ggplot themes to spice up your plots
for presentations, but some of these addtional features might be 
distractive and not be ideal for publications 
}
%\usage{
%
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%
%}
%\details{
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen
}
%\note{
%
%}
\seealso{
\code{\link{ggplot}}
}
\examples{
\dontrun{
data(ITS1, ITS2, meta)
data <- list(ITS1=ITS1, ITS2=ITS2)
# dissim.alleig.plot returns a ggplot2 object:
my.eig.plot <- dissim.alleig.plot(data)
my.eig.plot # view the plot

# update ggplot theme
require("grid") 
new_theme <-RAM.color()
my.eig.plot <- my.eig.plot + new_theme
my.eig.plot
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
