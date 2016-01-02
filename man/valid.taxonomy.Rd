\name{valid.taxonomy}
\alias{valid.taxonomy}
\alias{reformat.taxonomy}
\title{
Validate And Reformat The OTU Taxonomy Column
}
\description{
A properly formatted taxonomy column of an otu table is 
critical for \code{RAM} functions to run properly. The taxonomy 
column of an otu table is composed of taxonomic lineages for 
otuIDs. \code{RAM} accept 7 ranks, including kingdom, phylum, 
class, order, family, genus and species, sub ranks are not
supported. Taxa names at each rank should have prefix as "k__",
"p__", "c__", "o__", __", "g__", and "s__", each rank should be
separated by "; ", i.e. a semi colon and a white space, NOT 
just ";". This function will check the format of the taxonomy 
column of the input otu table and give suggetions that whether 
or not it needs to be reformatted using 
\code{\link{reformat.taxonomy}} of \code{RAM}. However, 
\code{RAM} does accept missing ranks in lineages, as long
as each rank is separated by "; " with proper prefix.
}
\usage{
valid.taxonomy(data)
reformat.taxonomy(data, input.ranks, sep="; ")
}
\arguments{
  \item{data}{
a list of otu tables, see RAM.input.formatting.
}
  \item{input.ranks}{
the ranks of the taxonomic lineages in the input otu tables.
}
  \item{sep}{
the separator between each taxonomic rank in the lineage.
}
}
\author{
Wen Chen.
}
\examples{
data(ITS1)
\dontrun{
# for demonstration purpose only
# the ITS1 dataset missing species rank, but it's ok
# the problematic taxonomy lineages are those missing proper 
# prefix at each rank
# see ?RAM.rank.formatting
valid.taxonomy(data=list(ITS1=ITS1))
input.ranks <- c("kingdom", "phylum", "class", "order", 
                 "family", "genus")
reform.ITS1 <- reformat.taxonomy(list(ITS1=ITS1),
                                 input.ranks=input.ranks,
                                 sep="; ")[[1]]
}
}
\seealso{
\code{\link{get.rank}}, \code{\link{tax.abund}}
}
