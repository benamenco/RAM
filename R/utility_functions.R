percent.classified <- function(data, rank) {
  valid.OTU(data)
  .valid.rank(rank)
  
  # number of rows in get.rank gives us the number of OTUs classified at that
  # rank; divide by total number of OTUs and multiply to get %
  return(100 * dim(get.rank(otu1=data, rank=rank))[1] / dim(data)[1])
}
