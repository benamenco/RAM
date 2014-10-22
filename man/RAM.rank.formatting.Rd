\name{RAM.rank.formatting}
\alias{RAM.rank.formatting}

\title{Taxonomic Rank Formatting}

\description{
In all \code{RAM} functions requiring the user to input a taxonomic 
rank, three different formats for this taxon are accepted. All of 
these formats are simple character vectors (strings), and are 
provided for readability and convenience. 
The user only needs to specify any single element from any of the 
formats below. The formats are as follows:
\describe{
\item{Format 1:}{"kingdom", "phylum", "class", 
                 "order", "family", "genus", "species"}
\item{Format 2:}{"k", "p", "c", "o", "f", "g", "s"}
\item{Format 3:}{"k__", "p__", "c__", "o__", 
                 "f__", "g__", "s__"}
}
}

\author{
Wen Chen and Joshua Simpson.
}

\seealso{
\code{\link{get.rank}}, \code{\link{tax.abund}}
}
