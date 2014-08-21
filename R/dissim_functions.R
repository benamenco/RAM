# These functions all deal with calculating measures related to dissimilarity matrices

dissim.clust <- function(data, dist.method="morisita", clust.method="average") {
  valid.OTU(data)
  otu.t <- transpose.OTU(data)
  
  # calculate the distances, then create the hclust
  dist <- vegdist(otu.t, method=dist.method)
  dist.clust <- hclust(dist, method=clust.method)
  
  return(dist.clust)
}

dissim.eig <- function(data, method="morisita") {
  valid.OTU(data)
  
  otu.t <- transpose.OTU(data)
  k.max <- dim(otu.t)[1] - 1
  
  dist <- vegdist(otu.t, method=method)
  # note: we suppress warnings because we will often have negative eigenvalues
  # due to numeric error (which raises a warning in pco) 
  dist.pco <- suppressWarnings(pco(dist, k=k.max))
  
  return(dist.pco$eig)
}

dissim.ord <- function(data, dist.method="morisita", k=NULL) {
  valid.OTU(data)
  
  k.max <- dim(data)[2] - 2
  
  if (is.null(k)) {
    k <- k.max
  }
  
  if (!is.numeric(k) || length(k) != 1L) {
    stop("k must be a numeric vector of length 1")
  }
  
  otu.t <- transpose.OTU(data)
  
  if (k <= k.max) {
      # calculate the ordination distances
      ord <- pco(vegdist(otu.t, method=dist.method), k=k)
      ord.dist <- vegdist(ord$points, method="euclidean")
      met.dist <- vegdist(otu.t, method=dist.method)
      return(list(ord.dist, met.dist))
  } else { 
    stop(paste("specified value for k was", k, 
               "which is larger than the maximum value", k.max)) 
  } 
}

dissim.GOF <- function(data, method="morisita") {
  valid.OTU(data)
  
  otu.t <- transpose.OTU(data)
  k.max <- dim(otu.t)[1] - 1
  
  met.dist <- vegdist(otu.t, method=method)
  
  GOF.vals <- numeric(k.max - 1)
  
  # calculate goodness of fit values for all possible values of k 
  # note: we suppress warnings because we will often have negative eigenvalues
  # due to numeric error (which raises a warning in cmdscale) 
  for (k in 2:k.max) {
    coord <- suppressWarnings(cmdscale(met.dist, k, eig=TRUE))
    GOF.vals[k] <- coord$GOF[2]
  }
  
  return(GOF.vals)
}

dissim.tree <- function(data, dist.method="morisita", clust.method="average") {
  valid.OTU(data)
  otu.t <- transpose.OTU(data)
  
  met.dist <- vegdist(otu.t, method=dist.method)
  # calculate the clustering and tree distances
  clust <- hclust(met.dist, method=clust.method)
  tree.dist <- cophenetic(clust)
  
  return(list(tree.dist, met.dist))
}

dissim.pvar <- function(data, method="morisita") {
  valid.OTU(data)
  
  otu.t <- transpose.OTU(data)
  k.max <- dim(otu.t)[1] - 1
  
  met.dist <- vegdist(otu.t, method=method)
  # note: we suppress warnings because we will often have negative eigenvalues
  # due to numeric error (which raises a warning in pco) 
  met.pcoa <- suppressWarnings(pco(met.dist, k=k.max))
  
  # get the number of positive eigenvalues
  axes <- dim(met.pcoa$points)[2]
  
  met.pcoa$eig / sum(met.pcoa$eig)
}
