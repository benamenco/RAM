\name{envis.NB}
\alias{envis.NB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize The Negative Binomial Model OF A Given Taxon OR OTUID
}
\description{
This function plot the negative binomial model for a given otuID or taxon
}
\usage{
envis.NB(NB.model="", tax.meta,  taxon="", 
         x="", num.col=NULL, group=NULL, group.order=NULL, 
         xlab=NULL, ylab=NULL, fill=NULL, facet=NULL, 
         file=NULL, ext=NULL, width=8, height=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NB.model}{
the negative binomial model. Can be obtained by using
\code{\link{assist.NB}}
}
  \item{tax.meta}{
the combined taxon/otuID and metadata. Can be obtained by using
\code{\link{assist.NB}}
}
  \item{taxon}{
the taxon or otuID. Can be obtained by using \code{\link{assist.NB}}
}
  \item{x}{
a metadata variable name for x axis.
}
  \item{num.col}{
optional. A metadata numerical variable that will be used as predictor.
}
  \item{group}{
optional. A metadata factor variable that will be used as predictor.
}
  \item{group.order}{
optional. The desired order for the \code{group}.
}
  \item{xlab}{
optional. X axis label.
}
  \item{ylab}{
optional. Y axis label.
}
  \item{fill}{
optional. Color for fill different categories.
}
  \item{facet}{
optional. Metadata category variables as faceting variables. 
}
  \item{file}{
optional. Filename that the plot to be saved to.  
}
  \item{ext}{
optional. Filename extension, type of image to be saved.  
}
  \item{width}{
an integer. Filter OTU table by counts.  
}
  \item{height}{
an integer. Filter OTU table by counts.  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function plot predicted taxon/otuID under the impact of metadata
variables.
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
data(ITS1, meta)

# filter otu table
its1 <- filter.OTU(data=list(ITS1=ITS1), percent=0.01)[[1]]
m <- meta[, c(2,3,5,7)]

\dontrun{
# test the model
nb <- assist.NB(its1, meta=m, rank="g", anov.fac="Harvestmethod", 
              taxon=rownames(its1)[1])

NB.model<-nb[[1]]
tax.meta<-nb[[2]]
taxon<-nb[[3]]

envis.NB(NB.model=NB.model, tax.meta=tax.meta,  taxon=taxon, 
                x="Ergosterol_ppm", num.col="Ergosterol_ppm", 
                group="Crop", group.order=NULL, 
                xlab="Ergosterol (ppm)", ylab=NULL, 
                fill="Harvestmethod", facet=c("City","Crop"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
