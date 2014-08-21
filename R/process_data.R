transpose.OTU <- function(data) {
  valid.OTU(data)
  # return the transpose (without the taxonomy column, which should be the last column)
  return(as.data.frame(t(data[, -dim(data)[2]])))
}

shared.OTU <- function(data) {
  valid.OTU(data)
  ###check for zero entries?
  
  otu <- transpose.OTU(data)
  otu.pa <- decostand(otu, "pa")

  # for readability
  x <- dim(otu.pa)[1]
  y <- dim(otu.pa)[2]
  
  num.otu.one.sample <- dim(otu.pa[ ,colSums(otu.pa) == 1])[2]
  num.otu.mult.sample <- dim(otu.pa[ ,colSums(otu.pa) > 1])[2] ### test this
  num.otu.all.sample <- dim(otu.pa[ ,colSums(otu.pa) == x])[2]
  
  per.otu.one.sample <- num.otu.one.sample / y
  per.otu.all.sample <- num.otu.all.sample / y
  
  num.seq.shared.otu <- sum(otu[ ,colnames(otu) %in% colnames(otu.pa[ ,colSums(otu.pa) == x])])
  per.seq.shared.otu <- num.seq.shared.otu / sum(otu)
  
  # get the OTU #s for all OTUs present in all samples
  otu.IDs <- colnames(otu.pa[ ,colSums(otu.pa) == x])
  # make a list of the form "OTU-taxonomic_information"
  tax <- paste(otu.IDs, "-", data[otu.IDs, ]$taxonomy, sep="")
  
  val <- list(num.otu.one.sample, num.otu.mult.sample, num.otu.all.sample,
           per.otu.one.sample, per.otu.all.sample, num.seq.shared.otu,
           per.seq.shared.otu, tax)
  
  # format numerical data?
  #val <- format(round(val[1:7], 2), nsmall=2)
  
names(val) <- c("#_of_OTUs_in_1_sample", "#_of_OTUs_in_>1_sample", 
                "#_of_OTUs_in_all_samples", "%_of_OTUs_in_one_sample",
                "%_of_OTUs_in_all_samples", "#_of_sequence_in_shared_OTUs",
                "%_of_sequence_in_shared_OTUs", "OTUs_in_all_samples")

  return(val)
}