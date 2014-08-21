#########################
# dissim.plot functions #
#########################

### A slightly ambitious task would be to abstract almost all the dissim.X.plot methods
### into one; the only (minor) difficulty is adjusting the plot settings for the different
### outputs

# small internal function to store the default distance methods
.get.dist.methods <- function() {
  return(c("morisita", "bray", "jaccard", "chao", "euclidean"))
}

# small internal function to store the default clustering methods
.get.clust.methods <- function() {
  return(c("average", "centroid", "complete", "mcquitty"))
}

.get.ord.methods <- function() {
  return(c("euclidean"))
}

dissim.clust.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL, 
                              clust.methods=NULL) {
  
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  given.dist.methods <- !is.null(dist.methods)
  given.clust.methods <- !is.null(clust.methods)
  

  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (!given.clust.methods) {
    clust.methods <- .get.clust.methods()
  }
  
  ncols <- length(dist.methods) * length(clust.methods)
  
  # if too many methods supplied, warn the user
  if (ncols >  10) {
    warning("you specified many distance and cluster methods; some may not fit on the plot. If none are displayed, try again with fewer.")
  }
  
  if (save) { 
    # the default size is too small; so adjust size manually
    .get.dev(file, "tiff", height = 5 * num.otus, width = 3 * ncols)
  }
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # determine the number of rows needed for plotting (one row for otu1, one
  # for otu2)
  par(mfcol=c(num.otus, ncols))
  
  # plot each distance/clustering method, for otu1 and (if given) otu2
  for (dist.met in dist.methods) {
    for (clust.met in clust.methods) {
      title <- paste(dist.met, "/", clust.met, sep="")
      
      # calculate the cluster data (also does data validation)
      dist.clust <- dissim.clust(otu1, dist.met, clust.met)
      # plot cluster data; hang=-1 aligns the labels
      plot(dist.clust, hang=-1, main=title)
      
      if (!single.otu) {
        # repeat for otu2
        dist.clust <- dissim.clust(otu2, dist.met, clust.met)
        
        plot(dist.clust, hang=-1, main=title)
      }
    }
  }
  
  if (save) { dev.off() }
  
  invisible()
}

dissim.eig.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL) {
  
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (save) { .get.dev(file, "tiff") }
  
  par(mfcol=c(num.otus, length(dist.methods)))
  
  for (met in dist.methods) {
    # calculate the eigenvalues (also does data validation)
    eigenvals <- dissim.eig(otu1, met)
    barplot(eigenvals, ylab="Eigenvalue", xlab="Axis Number", main=met)
    
    if (!single.otu) {
      eigenvals <- dissim.eig(otu2, met)
      barplot(eigenvals, ylab="Eigenvalue", xlab="Axis Number", main=met)
    }
  }
  
  if (save) { dev.off() }
  
  invisible()
}

dissim.ord.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL, 
                            k=NULL) {
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  if (single.otu) {
    k.max <- dim(otu1)[2] - 2
  } else {
    k.max <- min(dim(otu1)[2] - 2, dim(otu2)[2] - 2)
  }
  
  if (is.null(k)) {
    k <- k.max
  }
  
  if (!is.numeric(k) || length(k) != 1L) {
    stop("k must be a numeric vector of length 1")
  }
  
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  ncols <- length(dist.methods)
  
  # if too many methods supplied, warn the user
  if (ncols >  10) {
    warning("you specified many distance and ordination methods; some may not fit on the plot. If none are displayed, try again with fewer.")
  }
  
  if (save) { .get.dev(file, "tiff") }
  
  par(mfcol=c(num.otus, ncols))
  
  # plot the data, recall that output of dist.ord is a list where the first 
  # item is the ordination distances, the second is the given method distances
  for (dist.met in dist.methods) {
      title <- paste("dist-", dist.met, sep="")
      
      # calculate the distances (also does data validation)
      distances <- dissim.ord(otu1, dist.met, k)
      
      plot(distances[[1]], distances[[2]], 
           xlab="Ordination Distance", ylab="Method Distance", main=title)
      
      if (!single.otu) {
        distances <- dissim.ord(otu2, dist.met, k)
        
        plot(distances[[1]], distances[[2]], 
             xlab="Ordination Distance", ylab="Method Distance", main=title)
    }
  }
  
  if (save) { dev.off() }

  invisible()
}

dissim.GOF.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL) {
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (save) { .get.dev(file, "tiff") }
  
  par(mfcol=c(num.otus, length(dist.methods)))
  
  for (met in dist.methods) {
    # calculate the goodness of fit values (also does data validation)
    GOF.vals <- dissim.GOF(otu1, met)
    plot(GOF.vals, xlab="Dimensions", ylab="Goodness of Fit", main=met)
    
    if (!single.otu) {
      GOF.vals <- dissim.GOF(otu2, met)
      plot(GOF.vals, xlab="Dimensions", ylab="Goodness of Fit", main=met)
    }
  }
  
  if (save) { dev.off() }
  
  invisible()
}

dissim.tree.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL, 
                             clust.methods=NULL) {
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  given.dist.methods <- !is.null(dist.methods)
  given.clust.methods <- !is.null(clust.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (!given.clust.methods) {
    clust.methods <- .get.clust.methods()
  }
  
  ncols <- length(dist.methods) * length(clust.methods)
  
  # if too many methods supplied, warn the user
  if (ncols >  10) {
    warning("you specified many distance and cluster methods; some may not fit on the plot. If none are displayed, try again with fewer.")
  }
  
  if (save) { .get.dev(file, "tiff", width=3*ncols) }
  
  par(mfcol=c(num.otus, ncols))
  
  for (dist.met in dist.methods) {
    for (clust.met in clust.methods) {
      title = paste(clust.met, "/", dist.met, sep="")
      
      # calculate the distances (also does data validation)
      distances <- dissim.tree(otu1, dist.method=dist.met, clust.method=clust.met)
      # plot the data, recall that output of dist.tree is a list where the first 
      # item is the tree distances, the second is the given method distances
      plot(distances[[1]], distances[[2]], xlab="Tree Distances", 
           ylab="Method Distances", main=title)
      
      if (!single.otu) { # plot for otu2
        distances <- distances <- dissim.tree(otu2, dist.method=dist.met, 
                                              clust.method=clust.met)
        
        plot(distances[[1]], distances[[2]], xlab="Tree Distances", 
             ylab="Method Distances", main=title)
      }
    }
  }
  
  if (save) { dev.off() }
  
  invisible()
}

dissim.pvar.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL) {
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  given.dist.methods <- !is.null(dist.methods)
  
  # save par settings, restore when done plotting
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  if (save) { .get.dev(file, "tiff") }
  
  par(mfcol=c(num.otus, length(dist.methods)))
  
  for (met in dist.methods) {
    # calculate the percent variation (also does data validation)
    pct.var <- dissim.pvar(otu1, met)
    barplot(pct.var, xlab="Axis Number", ylab="% of Variation", main=met)
    
    if (!single.otu) {
      pct.var <- dissim.pvar(otu2, met)
      barplot(pct.var, xlab="Axis Number", ylab="% of Variation", main=met)
    }
  }
  
  if (save) { dev.off() }

  invisible()
}

dissim.alleig.plot <- function(otu1, otu2=NULL, file=NULL, dist.methods=NULL) {
  valid.OTU(otu1, otu2)
  
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  given.dist.methods <- !is.null(dist.methods)
  
  # if the user does not supply any methods, set the default ones
  if (!given.dist.methods) {
    dist.methods <- .get.dist.methods()
  }
  
  # parameter for pco calculations; can be at most n-1 where n is the # of samples
  k.max = dim(otu1)[1] - 1

  df.rows <- vector(mode="list", length=num.otus * length(dist.methods))
  
  index <- 1
  
  for (data in list(otu1, otu2)) {
    # stop if we don't have otu2 data
    if (is.null(data)) {
      break
    }
    
    # set region for rows
    ### can this be improved?
    if (index == 1) {
      region <- "otu1"
    } else {
      region <- "otu2"
    }
    
    for (met in dist.methods) {
      # get the eigenvalules
      eigenvals <- dissim.eig(data, met)
      # we want to save the fraction of the sum for each eigenvalue
      df.rows[[index]] <- c(met, region, eigenvals / sum(eigenvals))
      index <- index + 1
    }
  }
 
  # make our data frame; set stringsAsFactors must be false, otherwise our 
  # numeric data becomes factors (this way the numeric data becomes character
  # data which we can easily convert back)
  eigs <- as.data.frame(do.call(rbind, df.rows), stringsAsFactors=FALSE)
  
  # the number of columns with eigenvalue data
  eig.cols <- dim(eigs)[2]-2
  
  # convert our eigenvalues measures back to numeric data
  for (i in 3:eig.cols) {
    eigs[ ,i] <- as.numeric(eigs[ ,i])
  }
  
  # make the column names "Method", "Region", "1", ..., "n"
  colnames(eigs) <- c("Method", "Region", paste(1:eig.cols))
  
  eigs <- melt(eigs, id.vars=c("Method", "Region"))
  eigs <- rename(eigs, list(variable="eigenval", value="pct_of_sum"))
  eigs$pct_of_sum <- as.numeric(eigs$pct_of_sum)
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  p <- ggplot(eigs, aes_string(x="eigenval", y="pct_of_sum", group="Method", 
                               colour="Method")) +
       scale_y_continuous(labels = percent_format()) + 
       ylab("Fraction of Sum") +
       xlab("Axis Number") +
       facet_wrap(~Region) +
       geom_line(size=1) +
       geom_point(size=2)
  
  if (save) {
    file <- .ensure.filepath(file, "tiff")
    ggsave(filename=file)
  }
  
  p
}
