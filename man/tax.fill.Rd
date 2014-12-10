\name{tax.fill}
\alias{tax.fill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fill Missing Taxonomic Information
}
\description{
This function consumes an OTU table and returns a new OTU table 
where the taxonomic classifications which are unidentified, 
unclassified, incertae sedis, or simply missing, are replaced 
with a more descriptive entry.
}
\usage{
tax.fill(data, downstream = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the OTU table to be used.
}
  \item{downstream}{
logical. Should the replacement occur downstream, or upstream? 
(see Details)
}
}
\details{
If \code{downstream == TRUE}, the function will start at the 
kingdom level and work its way down. Whenever an invalid group is 
encountered (i.e. one of "unclassified", "unidentified", 
"incertae_sedis", or simply missing, ignoring capitalization), 
the last known 'good' group will be substituted in the form 
"\code{p__belongs_to_k_Fungi}." If \code{downstream == FALSE}, 
the function will begin at the species level and work up.

This example should help clarify: given the taxonomy "k__Fungi; 
p__unidentified; c__Tremellomycetes", the new taxonomy is as 
follows (recall that an OTU table is required as input, and will 
be returned as output; this example simply shows the effect on 
the taxonomy):
  \itemize{
    \item Downstream (Kingdom -> Species): "k__Fungi; 
    p__belongs_to_k_Fungi; c__Tremellomycetes; 
    o__belongs_to_c_Tremellomycetes; 
    f__belongs_to_c_Tremellomycetes; 
    g__belongs_to_c_Tremellomycetes;
    s__belongs_to_c_Tremellomycetes"
    
    \item Upstream (Species -> Kingdom): "k__Fungi;
    p__belongs_to_c_Tremellomycetes; 
    c__Tremellomycetes; 
    o__belongs_to_no_taxonomy; 
    f__belongs_to_no_taxonomy; 
    g__belongs_to_no_taxonomy; 
    s__belongs_to_no_taxonomy"
  }

Usually, \code{downstream = TRUE} will provide a more useful 
output, however if the species is often known for your data, but 
other ranks are unknown, \code{downstream = FALSE} will provide 
a more descriptive taxonomy.
}
\value{
Returns an OTU table as a data frame with the exact same numerical 
counts as \code{data}, but an updated \code{taxonomy} column.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Wen Chen and Joshua Simpson.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }
%% {....} ~

\seealso{
\code{\link{RAM.rank.formatting}}
}
\examples{
data(ITS1)

#\code{filter.OTU} returns a list
otu <- filter.OTU(data=list(ITS1=ITS1), percent=0.001)[[1]]

tax.fill(otu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
