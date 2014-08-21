get.rank <- function(otu1, otu2=NULL, rank=NULL) {
  single.otu <- is.null(otu2)
  single.rank <- !is.null(rank)
  valid.OTU(otu1, otu2)
  
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # if given both otu1 and otu2, call get.rank for both
  if (!single.otu) {
    
    output <- list()
    output$otu1 <- get.rank(otu1=otu1, rank=rank)
    # this looks ugly, but we are just calling get.rank with a single
    # OTU argument (which is named otu2)
    output$otu2 <- get.rank(otu1=otu2, rank=rank)
    
    if (single.rank) {
      names(output$otu1) <- .get.rank(.get.rank.ind(rank))
      names(output$otu2) <- .get.rank(.get.rank.ind(rank))
    } else {
      names(output$otu1) <- tax.classes
      names(output$otu2) <- tax.classes
    }
    
    return(output)
  }
  
  if (single.rank) {
    .valid.rank(rank)
    
    # select only the rows with the given taxonomic prefix in their taxonomy
    pat <- .get.rank.pat(rank)
    output <- otu1[grep(pat, otu1$taxonomy), ]
    
    # now only select the ones that are NOT on the blacklist
    remove.pat <- paste0(.blacklist(rank), "|no_taxonomy")
    
    output <- output[!grepl(remove.pat, output$taxonomy, ignore.case = TRUE), ]
    
    if (dim(output)[1] == 0) {
      warning(paste("no OTUs classified at the", .get.rank(.get.rank.ind(rank)),
                    "level."))
    }
    return(output)
    
  } else { # call get.rank for each taxonomic classification
    output <- vector(mode="list", length=length(tax.classes))
    
    for (i in tax.classes) {
      output[.get.rank.ind(i)] <- list(get.rank(otu1=otu1, rank=i))
    }
    
    names(output) <- tax.classes
    return(output)
  }
}

tax.split <- function(otu1, otu2=NULL, rank=NULL) {
  valid.OTU(otu1, otu2)
  single.otu <- is.null(otu2)
  single.rank <- !is.null(rank)
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # if given otu1 and otu2, call tax.split for both
  if (!single.otu) {
    
    output <- list()
    
    output$otu1 <- tax.split(otu1=otu1, rank=rank)
    output$otu2 <- tax.split(otu1=otu2, rank=rank)
    
    return(output)
  }
  
  if (single.rank) {
    # get the index for rank (also does data validation for rank)
    .valid.rank(rank)
    tax.ind <- .get.rank.ind(rank)
  }
  
  # split OTU table taxonomy information into individual columns
  otu.split <- concat.split(otu1, split.col="taxonomy", sep=";", drop=TRUE)
  
  # columns from 1 to max contain the numeric data, the other taxonomic
  max <- dim(otu1)[2] - 1
  
  # we need seven taxonomy columns; if one is missing (because nothing classified
  # at that level), add empty column
  
  # while we have less columns than we should...
  while (dim(otu.split)[2] < max + length(tax.classes)) {
    # ...add a column filled with empty strings
    otu.split <- cbind(otu.split, rep("", times=dim(otu.split)[1]))
  }
  
  # rename the columns
  names(otu.split)[-(1:max)] <- tax.classes
  # strip the 'formatting' bits of taxonomic info
  otu.split[ ,-(1:max)] <- gsub("k__|p__|c__|o__|f__|g__|s__|;", "", 
                                as.matrix(otu.split[ ,-(1:max)]))
  
  if (single.rank) {
    # add the single taxonomic column to the numeric data, return that data frame
    return(otu.split[ , names(otu.split) %in% c(names(otu1)[1:max], tax.classes[tax.ind])])
    
  } else {
    # set up list for output
    otu.taxa <- vector("list", length(tax.classes))
    names(otu.taxa) <- tax.classes 
    
    for (i in 1:length(tax.classes)) {
      # creating a list of data frames is surprisingly difficult in R;
      # if you do not use the list(df1, df2, ...) syntax, you get a list composed
      # of the first data frame in its entirety, followed be the individual columns
      # of the other data frames. Instead we wrap each data frame in a list itself 
      # before we add it, then we can access them with [[]]
      otu.taxa[i] <- list(otu.split[ , names(otu.split) %in% c(names(otu1)[1:max],
                                                               tax.classes[i])])
    }
    
    return(otu.taxa)
  }
}

tax.abund <- function(otu1, otu2=NULL, rank=NULL, drop.unclassified=FALSE, 
                      top=NULL, count=TRUE, mode="number") {
  
  single.otu <- is.null(otu2)
  valid.OTU(otu1, otu2)
  single.rank <- !is.null(rank)
  # data validation for top is done later in the function (when the dimensions of 
  # the taxonomy matrix are known)
  filter <- !is.null(top)
  
  if (!is.logical(drop.unclassified) || length(drop.unclassified) != 1L) {
    stop("argument 'drop.unclassified' must be a logical vector of length one.")
  }
  
  if (!is.logical(count) || length(count) != 1L) {
    stop("argument 'count' must be a logical vector of length one.")
  }
  
  if (!is.character(mode) || !any(mode %in% c("number", "percent"))) {
    stop("argument 'mode' must be one of 'number' or 'percent'.")
  }
  
  # if given otu1 and otu2, call tax.abund for both
  if (!single.otu) {
    
    drop <- drop.unclassified
    output <- list()
    
    output$otu1 <- tax.abund(otu1, rank=rank, drop.unclassified=drop, top=top,
                             count=count, mode=mode)
    output$otu2 <- tax.abund(otu2, rank=rank, drop.unclassified=drop, top=top,
                             count=count, mode=mode)
    
    return(output)
  }
  
  # get the OTU table in proper format
  if (single.rank) {
    .valid.rank(rank)
    tax.list <- list(tax.split(otu1, rank=rank))
  } else {
    tax.list <- tax.split(otu1)
  }
  
  for (i in seq(along=tax.list)) { 
    
    # update taxonomy label to "taxonomy"
    names(tax.list[[i]])[dim(tax.list[[i]])[2]] <- "taxonomy" 
    # aggregate the otu table by taxonomy rank names 
    tax.list[[i]] = aggregate(tax.list[[i]][ , -dim(tax.list[[i]])[2]],  
                              by = list(tax.list[[i]]$taxonomy), FUN = .sumcol) 
    # change row names to header provided by aggregate
    rownames(tax.list[[i]]) <- tax.list[[i]][ , 1] 
    # remove first column (that information is now in the row names)
    tax.list[[i]] <- tax.list[[i]][ , -1]
    # transpose table (the as.data.frame generates the V1 heading) 
    tax.list[[i]] <- as.data.frame(t(tax.list[[i]])) 
    
    # if count is false, return relative abundance
    if (!count) {
      tax.list[[i]] <- decostand(tax.list[[i]], method="total")
    }
    
    # order the table by column sums
    tax.list[[i]] <- tax.list[[i]][ , order(colSums(tax.list[[i]]), 
                                            decreasing = TRUE), drop=FALSE] 
    
    # remove all zero entries
    tax.list[[i]] <- tax.list[[i]][ , colSums(tax.list[[i]]) > 0, drop=FALSE]
    
    # rename the "V1" column
    if (!single.rank) {
      names(tax.list[[i]])[names(tax.list[[i]]) == "V1"] <- 
        paste("unclassified_below_", .get.rank(i - 1), sep="")
    } else {
      # if we are only processing one element, we cannot use the i index 
      # to get the correct rank
      names(tax.list[[i]])[names(tax.list[[i]]) == "V1"] <- 
        paste("unclassified_below_", .get.rank(.get.rank.ind(rank) - 1), sep="")
    }
    
    # remove unclassified columns
    if (drop.unclassified) {
      # this selects all columns NOT containing in the blacklist
      remove.pat <- paste0(.blacklist(.get.rank(i)), "|no_taxonomy")
      
      tax.list[[i]] <- tax.list[[i]][ , !grepl(remove.pat, names(tax.list[[i]]),
                                               ignore.case=TRUE), 
                                     drop=FALSE]
    }
    
    # keep only the 'top' most abundant groups, where top is user-given 
    if (filter) {
      tax.list[[i]] <- .select.top(tax.list[[i]], top, count, mode)
    }
  }
  
  if (single.rank) {
    return(tax.list[[1]])
    
  } else {
    return(tax.list)
  }
}

.select.top <- function(abund, top, count, mode) {
  if (!is.numeric(top) || length(top) != 1) {
    stop("argument 'top' must be a numeric vector of length one.")
  }
  
  if (!is.character(mode) || !any(mode %in% c("number", "percent"))) {
    stop("argument 'mode' must be one of 'number' or 'percent'.")
  }
  
  if (!is.numeric(top) || length(top) != 1L) {
    stop("argument 'top' must be a numeric vector of length one.")
  }
  
  if (top <= 0) {
    stop("argument 'top' must be greater than zero.")
  }
  
  # take top X samples when mode == "number"
  if (mode == "number") {
    num.groups <- dim(abund)[2]
    
    if (top > num.groups) {
      warning("argument 'top' is greater than the total number of groups, returning all groups.")
      top <- num.groups
    }
    
    abund <- abund[ ,1:top, drop=FALSE]
    
    
  # find all samples with RA > X% when mode == "percent"
  } else if (mode == "percent") {
    
    if (top > 100) {
      stop("argument 'top' must be in the range (0, 100].")
    }
    
    percent <- top / 100
    
    # we need the 'count' parameter to determine whether we need to standardize
    # the data ourselves or not
    if (count) {
      abund.stand <- decostand(abund, "total")
    } else {
      abund.stand <- abund
    }
    
    abund.max <- apply(abund.stand, MARGIN=2, FUN=max)
    # find samples w/ max relative read abundance < 'top'% and remove
    exclude <- names(which(abund.max < percent))
    
    # if all groups have been excluded
    if (length(exclude) == dim(abund)[2]) {
      stop("no taxon groups with relative abundance above 'top' percent.")
    }
    
    # select the samples above 'top'%
    abund <- abund[ ,-(which(names(abund.max) %in% exclude)), drop=FALSE]
  }
  
  abund
}

tax.fill <- function(data, downstream=TRUE) {
  
  valid.OTU(data)
  
  if (length(downstream) != 1L || !is.logical(downstream) || is.na(downstream)) {
    stop("'downstream' must be a character vector of length one.")
  }
  
  if (downstream) {
    range <- 1:7
  } else {
    range <- 7:1
  }
  
  taxonomy <- str_split(data$taxonomy, "; ")
  taxonomy.length <- length(taxonomy)
  
  # initialize the vector of classified names
  classified.names <- rep("no_taxonomy", times=taxonomy.length)
  # initialize list for fixed taxonomy
  new.taxonomy <- vector(mode="list", length=taxonomy.length)
  
  for (i in range) {
    # this selects the i-th element of every list and returns a character vector
    current.groups <- unlist(lapply(taxonomy, FUN="[", i))
    # any entries that are TRUE (i.e. match blacklist) need to be replaced
    # we add that pattern to match any group prefix with no name 
    # (since we split on semicolons the rank argument in .blacklist won't work)
    blacklist <- paste0(.blacklist(), "|[kpcofgs]__$")
    
    replace <- is.na(current.groups) | grepl(blacklist, current.groups,
                                             ignore.case = TRUE)
    
    num.groups <- length(current.groups)
    for (j in 1:num.groups) {
      
      # if we need to replace the classification, use last good name
      if (replace[j]) {
        
        new.name <- paste0(.get.rank.pat(.get.rank(i)), "belongs_to_",
                           gsub("__", "_", classified.names[j]))
        
        new.taxonomy[[j]][i] <- new.name
      } else { # if the classification is good, store it
        classified.names[j] <- current.groups[j]
        new.taxonomy[[j]][i] <- classified.names[j]
      }
    }
  }
  
  new.taxonomy <- lapply(new.taxonomy, FUN=paste, collapse="; ")
  
  # return a new data frame with the updated taxonomy
  # (do not mutate the old data frame)
  data.frame(data[ ,-dim(data)[2], drop=FALSE], taxonomy=unlist(new.taxonomy),
             stringsAsFactors = FALSE)
}