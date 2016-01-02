\name{group.diversity}
\alias{group.diversity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boxplot To Compare Diversity Indices Among Groups
}
\description{
This function first use \code{\link{OTU.diversity}} to calculate
the diversity indices for each sample and then do a boxplot to
compare the selected indices among different groups.
}
\usage{
group.diversity(data, meta, factors="", indices="",
                diversity.info=FALSE,
                x.axis=NULL, compare=NULL,
                facet=NULL, facet.y=TRUE, facet.x.cex=NULL,
                facet.y.cex=NULL, scale.free=NULL,
                xlab=NULL, ylab=NULL,
                legend.title=NULL, legend.labels=NULL,
                file=NULL, ext=NULL, width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a list, containing otu tables. 
See also \code{\link{RAM.input.formatting}}
}
  \item{meta}{
the metadata table to be used (must have same samples as
\code{data}.
}
  \item{factors}{
a character vector. Must be variables in the metadata
}
  \item{indices}{
a character vector. Must be one or more of the following: 
"spec", "sim", "invsim", "shan", "sim_even", "shan_even", 
"sim_trudiv", "shan_trudiv", "chao", "ACE". See also 
\link[RAM]{OTU.diversity}, \code{\link{true.diversity}}, 
\code{\link{evenness}}, and \code{\link{diversity}}.
}
  \item{diversity.info}{
logical. Whether the diversity indices have calculated and 
included in the metadata table. The diversity indices should be 
processed by \code{\link{OTU.diversity}} for the same otu tables
and metadata table.
}
  \item{x.axis}{
optional. If \code{NULL}, will use the first variable in
\code{factors}; otherwise, must be one factor in the metadata or
'SampleID'
}
  \item{compare}{
optional. If \code{NULL}, will use the first variable in
\code{factors}; otherwise, must be one factor in the metadata
}
  \item{facet}{
optional. If provided, must be one factor in the metadata or 
'SampleID'
}
  \item{facet.y}{
logical, whether the \code{facet} being used as strip text of y
axis or x axis.
}
  \item{facet.x.cex}{
optional, an integer, the font size of the \code{stip.text.x} in
\code{ggplot}
}
  \item{facet.y.cex}{
optional, an integer, the font size of the \code{strip.text.y} 
in \code{ggplot}.
}
  \item{scale.free}{
optional. Whether use free scale for y axis.
}
  \item{xlab}{
optional. If not provided, the \code{x.axis} will be used as the
title of the x axis, otherwise, will use the provided string.
}
  \item{ylab}{
optional. If not provided, "value" will be used as the
title of the y axis, otherwise, will use the provided string.
}
  \item{legend.title}{
optional. If not provided, \code{compare} will be used as the
title of the legend, otherwise, will use the provided string.
}
  \item{legend.labels}{
optional. If not provided, will use the levels of \code{compare}
for the legends, otherwise, will use the provided vector of 
strings. The length of the provided vector of strings must 
equals to the levels of \code{compare}.
}
  \item{file}{
the filename to save the plot.
}
  \item{ext}{
the extention (file type) of the plot to saved.
}
  \item{width}{
the width of the plot to be saved.
}
  \item{height}{
the heigth of the plot to be saved.
}
}
%\details{
%%  Maybe explain full/restricted models
%}
\value{
This function returns a boxplot to compared selected diversity
indices among different groups.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen.
}
%\note{
%%  ~~further notes~~
%}
%% ~Warning with \section{Warning }{....} ~
\seealso{
\code{\link{OTU.diversity}}, \code{\link{true.diversity}},
\code{\link{evenness}} and \code{\link{diversity}}
}
\examples{
data(ITS1, ITS2, meta)
\dontrun{
RAM.theme<-RAM.color()
group.diversity(data=list(ITS1=ITS1, ITS2=ITS2), meta=meta,
                factors=c("Crop", "City"),
                indices=c("sim_trudiv", "shan_trudiv"),
                x.axis="Crop", compare="Harvestmethod",
                facet="City", facet.y=FALSE) + RAM.theme
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
