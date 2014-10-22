correlation<-function(data=NULL, is.OTU=TRUE, meta=NULL, rank="g", sel=NULL, 
                      sel.OTU=TRUE, data.trans=NULL, method="pearson", 
                      main=NULL, file=NULL, width=8, height=8) {

  save <- !is.null(file)
  if (save) { tiff(filename = file, compression="lzw", width=width, 
                   height=height, res=1000, units="in") }

  # make sure either data or metadata was provided
  if ( is.null(data) & is.null(meta) ) {
     stop("Error: provide at least one of the following: 
            OTU table, taxonomy abundance matrix, or metadata")
  } 

  #get rank names
  rank.name<-.get.rank.name(rank, plural=TRUE)
  rank_name<- .get.rank.name(rank)
  rank_pat<-.get.rank.pat(rank)

  if ( !is.null(meta) ) {
      # metadata of numeric variables only
      col <- sapply(meta, is.numeric)
      meta.sel <- meta[, col]
      if ( ncol(meta) == 0 ) {
          warning("no columns in metadata is numeric, will only plot data")
          meta.sel <- data.frame(matrix(vector(), nrow(meta.sel), 0))
      } else {
          meta.sel <- meta.sel
    }
  } else {
      if ( is.OTU ) {
          meta.sel <- data.frame(matrix(vector(), ncol(data)-1, 0))
      } else {
          meta.sel <- data.frame(matrix(vector(), nrow(data), 0))
      }
  }
 
  # validate data and metadata
  if ( !is.null(data) & !is.null(meta) ) {
      if ( is.OTU ) {
          valid.OTU(data)
          .valid.meta(otu1=data, meta=meta)         
      } else {
          data <- data[match(rownames(meta), rownames(data)), ]
          if ( !identical(rownames(data), rownames(meta)) ) {
              stop("Error: not same subjects in data and metadata")
          } 
      }
  } else if ( !is.null(data) & is.null(meta) ) {
      data <- data
  } else {
      data <- data
  }
  
  # require count data
  if ( !is.null(data) ) {
    # check whether integer data being provided
     if ( is.OTU ) {
          int <- unique(sapply(data[, -ncol(data)], is.integer))
     } else {
          int <- unique(sapply(data, is.integer))
     }
     if ( !isTRUE(int) & !is.null(data.trans) ) {
          warning("data contain non-integer records, are you sure to transform the data?")
      }
  }

  if ( !is.null(data) ) { 
      # plot OTU
      if ( sel.OTU ) {
          if ( !is.OTU ) {
              stop("Error: provided outIDs but data is not an OTU table")
          } else {
              if ( is.null(data.trans) ) {
                  data.sel <- data
              } else {
                  data.sel <- cbind(decostand(data[, -ncol(data)], 
                           MARGIN=2, data.trans), data$taxonomy)
                  names(data.sel)[ncol(data.sel)] <- "taxonomy"
              }
           }
           if ( !is.null(sel) ) {
               sel.otu<- which(rownames(data.sel) %in% sel)
               if ( length(sel.otu) != 0 ) {
                   data.sel <- data.sel[sel.otu, ]
               } else {
                   warning("did not find provided otuIDs in data,
                         use all data")
                   data.sel <- data.sel
               }
            } else {
               data.sel <- data.sel
            }
           
            data.tax <- LCA.OTU(data.sel, strip.format=FALSE, drop=TRUE)
            rownames(data.tax) <- paste( rownames(data.tax), 
                                    data.tax$LCA, sep="_")
            data.tax <- data.tax[, -ncol(data.tax)]
            data.tax <- as.data.frame(t(data.tax))
        }
  
      # plot taxa at given rank    
      if ( !sel.OTU ) {
          if ( is.OTU ) {
              data.tax <- tax.abund(data, count=TRUE, rank=rank, 
                                  drop.unclassified=FALSE)
          } else {
              data.tax <- data
          }
          # transform data
          if ( !is.null(data.trans) ) {
              data.tax <- decostand(data.tax, data.trans)
          } else {
              data.tax <- data.tax
          }
          # remove unclassified taxa
          remove.pat <- gsub(.get.rank.pat(rank), "", 
                       paste0(.blacklist(.get.rank.pat(rank)), "|no_taxonomy"))
          data.tax <- data.tax[, !grepl(remove.pat, names(data.tax), 
                        ignore.case=TRUE), drop=FALSE]
          # select only taxa in sel
          if ( !is.null(sel) ) {
              sel.tax <- which(colnames(data.tax) %in% sel)
              if ( length(sel.tax) != 0 ) {
                  data.tax <- data.tax[ , sel.tax]
              } else {
                  warning("No taxa in the provided list was found in data, 
                           will plot all taxa")
                  data.tax <- data.tax
              }
          } else {
              data.tax <- data.tax
          }
      } 
  } else {
      data.tax <- data.frame(matrix(vector(), nrow(meta), 0))
  }
  
  if (require("lattice")) {
      lattice::levelplot
      lattice::do.breaks
  } else {
    stop("package 'lattice' is required to use this function")
  #  source("http://bioconductor.org/biocLite.R"); biocLite("Heatplus")
  } 
  # combine data
  dat <- cbind(data.tax, meta.sel)

  if ( ncol(dat) == 0 ) {
      stop("Error: nothing to plot")
  } else if ( ncol(dat) > 200 ) {
      stop("Too many variables to be plotted, please use sel 
            option to reduce the number of variables")
  } else {
      
      # stats::cor
      # method <- c("pearson", "kendall", "spearman")
      cor.mat<-cor(dat, method=method, use = "pairwise.complete.obs")
      x.scale <- list(cex=0.5, alternating=1, col='black') 
      y.scale <- list(cex=0.5, alternating=1, col='black') 
  
      # lattice::levelplot; lattice::do.breaks
      bp1 <- lattice::levelplot(cor.mat,xlab=NULL,ylab=NULL,cex=0.08,  
                   at=lattice::do.breaks(c(-1.01,1.01),101),
                   scales=list(x=list(rot=90)),
                   colorkey=list(space="top"), 
                   col.regions=colorRampPalette(c("red","white","blue")), 
                   main=main,cex=1)
      print(bp1)
      
      if (save) { dev.off() }
              invisible()
  }
 
}

reset.META <- function(meta, factor=NULL, numeric=NULL, date=NULL) {

  if( !is.null(factor) ) {
    for ( i in factor ) {
      if ( ! (i %in% names(meta) ) ) {
        warning(paste (i," is not in in metadata!", sep=""))
      } else {
        meta[[i]] <- as.factor(meta[[i]])
      }
    }
    meta<-meta
  } else {
    meta<-meta
  }

  if( !is.null(numeric) ) {
    for ( i in numeric ) {
      if ( ! ( i %in% names(meta) ) ) {
         warning(paste (i," is not in in metadata!", sep=""))
      } else if ( !is.numeric(meta[[i]]) ) {
         warning(paste(i, " is not a numeric data type, but will
                       be converted to numeric as required!", sep=""))
         meta[[i]] <- as.numeric(as.factor(meta[[i]]))
      } else {
        meta[[i]] <- as.numeric(meta[[i]])
      }
    }
    meta<-meta
  } else {
    meta<-meta
  }

  if( !is.null(date) ) {
    for ( i in date ) {
      if( ! ( i %in% names(meta) ) ) {
        warning(paste (i," is not in in metadata!", sep=""))
      } else {
        meta[[i]] <- as.Date(meta[[i]])
      }
    }
    meta<-meta
  } else { 
    meta<-meta
  }
  
  return(meta)
}


