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
  
  if(length(which(colSums(otu.pa) == x))==0) {
    stop("NO OTUs are shared by all samples")
  } else {

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
}

filter.OTU <- function(otu1, otu2=NULL, percent=NULL, number=NULL) {

  if(is.null(otu2)) {
  	num.otus<-1
  } else {
	num.otus<-2
  }
	
  otu_sel <- list()
  index <- 1
  for (otu in list(otu1, otu2)){
  	# exit if otu2 is null
  	if (is.null(otu)) { break }
		
	valid.OTU(otu)

	# both filter are NULL
	if (is.null(number) & is.null(percent)) {

		warning("No filter was provided, all otus will be used. Consider filtering otu by total number of sequence (e.g. 50) or maximum relative abundance (e.g. 0.01)")
		sel<-1:nrow(otu)

	} else if (!is.null(number) & !is.null(percent)) {

		stop("Error: can only provide one type of filter by total number of sequence (e.g. 50) or maximum relative abundance (e.g. 0.01)")

	} else if (!is.null(percent)) {

		# select OTUs more than percent in at least on sample (otu columns)
		otu.p <- cbind(decostand(otu[, -ncol(otu)], MARGIN=2, "total"), otu$taxonomy)
		names(otu.p)[ncol(otu.p)]<-"taxonomy"
		sel <- which(apply(otu.p[, -ncol(otu.p)], MARGIN=1, FUN=max) > percent)
	} else if (!is.null(number)) {

		# select OTUs more than number sequences in total
		sel<-which(rowSums(otu[, -ncol(otu)]) > number)
	} else {
		
	} 

	otu_sel[[index]] <- sel
	index <- index + 1
  } 
	
  if (num.otus == 1) {
	if ( length(otu_sel[[1]]) == 0) {
		warning("no otus met the filter requirement, original OTU returned")
		return(otu1)
	} else {
		print(paste(length(otu_sel[[1]]), " otus met the filter requirment", sep=""))
		return(list(otu1[otu_sel[[1]], ], otu_sel[[1]]))
	}
  } 
  
  if (num.otus == 2) {
	otu.sel = list()

  	if ( length(otu_sel[[1]]) == 0 ) {
		warning("no otus in OTU1 met the filter requirement, original OTU returned")
		otu.sel[[1]] <- otu1
		otu_sel[[1]] <- ""
	} 
	if ( length(otu_sel[[2]]) == 0 ) {
		warning("no otus in OTU2 met the filter requirement, original OTU returned")
		otu.sel[[2]] <- otu2
		otu_sel[[2]] <- ""
	} 
	if ( length(otu_sel[[1]]) != 0 & length(otu_sel[[2]]) != 0 ) {
		print(paste(length(otu_sel[[1]]), " otus in OTU1 & ", length(otu_sel[[2]]), " otus in OTU2 met the filter requirment.", sep=""))
		otu.sel[[1]] <- otu1[otu_sel[[1]], ]
		otu.sel[[2]] <- otu2[otu_sel[[2]], ]
	}
	return(list(otu.sel[[1]], otu.sel[[2]], otu_sel[[1]], otu_sel[[2]]))
  }
 }


OTU.subsets<-function(otu1=otu1, otu2=NULL, meta=meta, factor="", value="", and=TRUE){

  if (is.null(otu2)) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  # meta factor and value for selection
  factors <- parse(text = paste("meta$",factor,sep=""))
  values <- paste(value, collapse="|")

  if( length(factor) > 0 ) {
      if (!isTRUE(and)) {
          # and=FALSE; any subjects belong to the value should be kept
          meta_sel <- list()
          for ( i in 1:length(factor) ) {
              meta_sel[[i]] <- meta[grepl(values, eval(factors[i])), ]
          }
          meta_sel <- meta_sel[lapply(meta_sel,length)>0]
          meta_sub <- unique(do.call("rbind", meta_sel))

      } else {
          # and=TRUE: subjects should belong to all values to be kept
          arg <- vector()
          for ( i in 1:length(factor) ) {
              new <- paste("grepl( \"", values, "\", ", "meta$", factor[i], " )", sep="")
              arg <- c(arg, new)
          } 
          arg.all <- paste(arg, collapse = " & ")
          sel <- eval(parse(text = arg.all))        
          meta_sub <- meta[sel, ]
          
    }
  }

  if ( nrow(meta_sub) == 0 ) {
      stop("Error: no subject met the filter requirement")
  } else {
      for (i in factor) {
          if (is.factor(meta[[i]])) {
              # drop levels
              meta_sub[[i]] <- factor(meta_sub[[i]])
          }
      }

     # match samples among OTUs and metadata
      result=list()
      for (i in 1:num.otus) {
          otu=list(otu1,otu2)[[i]]
          otu_sub = otu[, match(c(rownames(meta_sub),"taxonomy"), colnames(otu))]  
          ### remove empty rows and columns and order the data.frame ###
          otu_sub = otu_sub[rowSums(otu_sub[,-ncol(otu_sub)])>0, colSums(otu_sub[,-ncol(otu_sub)])>0]
          otu_sub = otu_sub[order(rowSums(otu_sub[,-ncol(otu_sub)]),decreasing=TRUE),]  
          ### Verify if the data.frame is non empty ###
          if(nrow(otu_sub) == 0){
              return ("Length of factor data is 0, please check if there were errors in your factors and factor.values")
          }
          if (isTRUE(all(rownames(meta_sub)==colnames(otu_sub[,-ncol(otu_sub)])))) {
              result[[i]]<-otu_sub
          } else {
              stop("meta and `otu` do not have same samples")      
          }
    
       }
    if ( num.otus == 1 ) {
        return(list(result[[1]], meta_sub))
    } else {
        return(list(result[[1]], result[[2]], meta_sub))
    }
  }
}
