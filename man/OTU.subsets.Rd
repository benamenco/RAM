\name{OTU.subsets}
\alias{OTU.subsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset OTU  
}
\description{
This function subset OTUs and metadata based on user defined 
values of metadata variables.   
}
\usage{
OTU.subsets(otu1, otu2=NULL, meta, factor, value, and=TRUE)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu1}{
the first OTU table to be analyzed.  
}
  \item{otu2}{
optional. the second OTU table to be used.  
}
  \item{meta}{
the metadata for subset.  
}
  \item{factor}{
a vector containing metadata variables.     
}
  \item{value}{
a vector containing values of interest in metadata variables.    
}
  \item{and}{
logical. Determine whether all conditions needs to be met or not.    
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a list of length 2 or 3.  
\itemize{
\item If \code{otu2} is given: a list of length 3, containing the 
      following members: 1) subset of otu1 otu1; 2) subset of otu2; 
      3) subset of metadata.   
\item If \code{otu2} is not given: output is a list of length 2, 
      containing the following members: 1) subset of otu1 otu1; 
      2) subset of metadata.    
  }
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
    data(ITS1, ITS2, meta)
    names(meta)
    factor <- c("City", "Harvestmethod")
    value <- c("City1", "Method1")

    sub <- OTU.subsets(otu1=ITS1, otu2=ITS2, meta=meta, 
                           factor=factor, value=value, and=TRUE)
    length(sub)
    dim(sub[[1]])
    dim(sub[[2]])
    dim(sub[[3]])  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
