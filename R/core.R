core.OTU<-function(otu, meta, rank="g", meta.factor="", percent=1){

  # valid otu & meta
  valid.OTU(otu)
  #.valid.meta(otu1=otu, meta=meta)
  samp <- colnames(otu)[1:(ncol(otu)-1)]
  meta <- meta[match(samp, rownames(meta)), ]
  if ( !identical(rownames(meta), samp) ) {
    stop("Error: metadata and otu table do not have same samples")
  }

  # check meta category
 if( !length(meta.factor) == 1 ) {
    stop("Error: please provide one factor variable in metadata")
  } else if ( ! (meta.factor %in% names(meta)) ) {
    stop("Error: meta.factor is not in metadata table")
  } else if ( ! is.factor(meta[[meta.factor]]) ) {
    stop("Error: need qualitative variable of metadata")
  } else {
    fac.levels<-levels(factor(meta[[meta.factor]]))
  }

  # get rank name & pattern
  rank.name<-.get.rank.name(rank, plural=TRUE)
  rank_name<- .get.rank.name(rank)
  rank_pat<-.get.rank.pat(rank)
  
  # tax.split
  otu<-tax.split(otu1=otu, rank=rank)
  otu<-otu[otu[[rank_name]] !="", ]

  ###check for zero entries?
  otu.fac.list=list()
  otu.fac.id = list()
  otu.fac.rank = list()
  for (i in 1:length(fac.levels)) {
    # each fac levels
    meta.fac <- meta[grep(fac.levels[i], meta[[meta.factor]]), ]
    otu.fac <- cbind(otu[,match(rownames(meta.fac), colnames(otu))], 
                     otu[[rank_name]])
    colnames(otu.fac)[ncol(otu.fac)]<-rank_name
    otu.fac <- otu.fac[rowSums(otu.fac[, -ncol(otu.fac)])>0, 
                      colSums(otu.fac[, -ncol(otu.fac)])>0]

    # presence/absence
    otu.fac.pa <- decostand(otu.fac[, -ncol(otu.fac)], "pa")
    
    # select only otus present in more than percent of samples in each category
    sel <- which(rowSums(otu.fac.pa) >= round(percent*ncol(otu.fac.pa)))
    if ( length(sel) != 0) {
        otu.fac.ids <- as.character(rownames(otu.fac)[sel])
        otu.fac.id[[i]] <- as.character(rownames(otu.fac)[sel])
        otu.fac.ranks <- unique(as.character(factor(otu.fac[[rank_name]][sel])))
        otu.fac.rank[[i]] <- unique(as.character(factor(otu.fac[[rank_name]][sel])))

        otu.fac.ids.num <- length(otu.fac.ids)
        otu.fac.ranks.num <- length(otu.fac.ranks)
        otu.fac.ids.perc <-  100*sum(otu.fac[c(otu.fac.ids), -ncol(otu.fac)])/sum(otu.fac[, -ncol(otu.fac)])

    } else {
        otu.fac.ids <- ""
        otu.fac.id[[i]] <- ""
        otu.fac.ranks <- ""
        otu.fac.rank[[i]] <- ""
        otu.fac.ids.num <- 0
        otu.fac.ranks.num <- 0
        otu.fac.ids.perc <- 0
    } 
    # summarize
    description <- paste(otu.fac.ids.num, "_otus_found_in_", percent*100, "%_", 
                   fac.levels[i], "_samples; ","assigned_to_", otu.fac.ranks.num, "_", 
                   rank.name, "; ", otu.fac.ids.perc, "%_of_total_sequences_in_", 
                    fac.levels[i], sep="")
    otu.fac.list[[i]]<- list(summary=description, otuID=otu.fac.ids, taxa=otu.fac.ranks)
    names(otu.fac.list)[[i]] <- fac.levels[i]
  }

  # find otus across all categories        
  otu.fac.list.num <- length(otu.fac.list)

  otu.fac.id[[otu.fac.list.num+1]] <- as.character(factor(Reduce(intersect, otu.fac.id)))
  otu.fac.rank[[otu.fac.list.num+1]] <- as.character(factor(Reduce(intersect, otu.fac.rank)))

  # otu ids in all samples
  # percent of sequences in core otus across all categories
  
  if ( isTRUE(length(otu.fac.id[[otu.fac.list.num+1]]) == 0) ||
        unique(otu.fac.id[[otu.fac.list.num+1]]) == ""  ) {
    otu.fac.id.all.num <- 0
    otu.fac.rank.all.num <- 0
    otu.fac.perc <- 0
  } else {
    otu.fac.id.all.num <- length(unique(otu.fac.id[[otu.fac.list.num+1]]))
    otu.fac.rank.all.num <- length(unique(otu.fac.rank[[otu.fac.list.num+1]]))
    otu.sel <- otu[rownames(otu) %in% otu.fac.id[[otu.fac.list.num+1]], ]
    otu.fac.perc <- 100*sum(otu.sel[, -ncol(otu.sel)])/sum(otu[, -dim(otu)[2]])
  } 

  description <- paste(otu.fac.id.all.num, "_otus_present_in_core_OTU_of_each_", 
                    meta.factor, "; ", "assigned_to_", otu.fac.rank.all.num, "_", 
                     rank.name, "; ", otu.fac.perc, "%_of_total_sequences",sep="")

  otu.fac.list[[otu.fac.list.num+1]] <- list(summary=description, 
                                             otuID=otu.fac.id[[otu.fac.list.num+1]], 
                                              taxa=otu.fac.rank[[otu.fac.list.num+1]])

  names(otu.fac.list)[[otu.fac.list.num+1]] <- meta.factor
    
  return(otu.fac.list)

}


core.Taxa<-function(data, is.OTU=FALSE, rank="g", meta, meta.factor="", percent=1){

  if( !length(meta.factor)==1 ) {
    stop("Error: please provide one factor variable in metadata")
  } else {
    fac.levels<-levels(factor(meta[[meta.factor]]))
  }

  #get rank names
  rank.name<-.get.rank.name(rank, plural=TRUE)
  rank_name<- .get.rank.name(rank)
  rank_pat<-.get.rank.pat(rank)

  if ( is.OTU ) {
    valid.OTU(data)
    .valid.meta(otu1=data, meta=meta)
        
    # tax.abund
    tax.all <- tax.abund(data, rank=rank, drop.unclassified=FALSE, count=TRUE)
  } else {
    data<-data[match(rownames(meta), rownames(data)), ]
    if (identical(rownames(data), rownames(meta))) {
        tax.all <- data
    } else {
        stop("Error: metadata and data don't have same subjects")
    }
 }

  # remove unclassified columns
  # this selects all columns NOT containing in the blacklist
  remove.pat <- gsub(rank_pat, "", paste0(.blacklist(rank_name),
                    "|no_taxonomy"))
  tax <- tax.all[ , !grepl(remove.pat, names(tax.all), 
                  ignore.case=TRUE), drop=FALSE]

  ###check for zero entries?
  tax.fac.list=list()
  tax.fac.id = list()

  for ( i in 1:length(fac.levels) ) {

    # each fac levels
    meta.fac <- meta[grep(fac.levels[i], meta[[meta.factor]]), ]
    tax.fac <- tax[match(rownames(meta.fac), rownames(tax)),]
    tax.fac <- tax.fac[rowSums(tax.fac)>0, colSums(tax.fac)>0]
    # presence or absence
    tax.fac.pa <- decostand(tax.fac, "pa")

    # select taxa present in more than percent of samples in each category
    sel <- which(colSums(tax.fac.pa) >= round(percent*nrow(tax.fac.pa)))

    if (length(sel) != 0) {
        tax.fac.ids <- as.character(colnames(tax.fac)[sel])
        tax.fac.id[[i]] <- as.character(colnames(tax.fac)[sel])
        tax.fac.ids.num <- length(tax.fac.ids)
        tax.fac.ids.perc <-  100*sum(tax.fac[, c(tax.fac.ids)])/sum(tax.fac)

    } else {
        tax.fac.ids <- ""
        tax.fac.id[[i]] <- ""
        tax.fac.ids.num <- 0
        tax.fac.ids.perc <- 0
    } 
    # summarize
    description <- paste(tax.fac.ids.num, "_", rank.name, "_found_in_", 
                              percent*100, "%_", fac.levels[i], "_samples; ", 
                              tax.fac.ids.perc, "%_of_total_sequences_in_", 
                              fac.levels[i], sep="")
    tax.fac.list[[i]]<- list(summary=description, taxa=tax.fac.ids)
    names(tax.fac.list)[[i]] <- fac.levels[i]
  }

  # find taxa across all categories        
  tax.fac.list.num <- length(tax.fac.list)
  tax.fac.id[[tax.fac.list.num+1]] <- as.character(factor(
                                         Reduce(intersect, tax.fac.id)))

  # taxa in all samples
  # percent of sequences in core taxa across all categories
  if( unique(tax.fac.id[[tax.fac.list.num+1]]) == "" || 
       length(unique(tax.fac.id[[tax.fac.list.num+1]])) == 0 ) {
    tax.fac.id.all.num <- 0
    tax.fac.perc <- 0
    tax.fac.id.all <- ""
  } else {
    tax.fac.id.all.num <- length(tax.fac.id[[tax.fac.list.num+1]])
    tax.fac.id.all <- tax.fac.id[[tax.fac.list.num+1]]
    tax.sel <- tax[, which(colnames(tax) %in% tax.fac.id.all)]
    tax.fac.perc <- 100*sum(tax.sel)/sum(tax)
  } 

  description <- paste(tax.fac.id.all.num, "_", rank.name, "_found_in_", 
                              percent*100, "%_of_samples_at_each_", meta.factor,
                              "; ", tax.fac.perc, "%_of_total_sequences", 
                              sep="")
  tax.fac.list[[tax.fac.list.num+1]] <- list(summary=description, taxa=tax.fac.id.all)
  names(tax.fac.list)[[tax.fac.list.num+1]] <- meta.factor
    
  return(tax.fac.list)

}


group.OTU <- function(otu, rank="g", otuIDs="", meta, meta.factor="", 
                 boxplot=TRUE, log=FALSE, main=NULL, file=NULL, ext=NULL, 
                 height=8, width=16) {
  # validate inputs
  valid.OTU(otu)
  .valid.meta(otu1=otu, meta=meta)
  rank.name<-.get.rank.name(rank, plural=TRUE)
  rank_name<- .get.rank.name(rank)
  rank_pat<-.get.rank.pat(rank)

    if(!length(meta.factor)==1) {
        stop("Error: please provide one factor variable in metadata")
    } else {
        fac.levels<-levels(factor(meta[[meta.factor]]))
    }

  save <- !is.null(file)
  if (save) { tiff(filename = file, compression="lzw", width=width, 
             height=height, res=1000, units="in") }  
  
  # get the groups
    if(!length(meta.factor)==1) {
        stop ("Error: please provide one factor")
    }
    
    if(length(otuIDs)==0) {
        stop ("Error: please provide list of otuIDs to plot")
    }

    # tax.split
    otu<-tax.split(otu1=otu, rank=rank)
    otu.p <- cbind(decostand(otu[, -ncol(otu)], "total", MARGIN=2), otu[[rank_name]])
    colnames(otu.p)[ncol(otu.p)] <- rank_name
    otu.p<-otu.p[otu.p[[rank_name]] !="", ]

    if( length(which(rownames(otu.p) %in% otuIDs))==0) {
        stop ("Error: no selected otuIDs present in the otu table")
    }
    
    otuIDs<-unique(otuIDs)
    otu.p.sel<-otu.p[which(rownames(otu.p) %in% c(otuIDs)), ]
    otu.sel<-otu[which(rownames(otu) %in% c(otuIDs)), ]

    otus.ex<-vector()
    for (i in otuIDs) {
        if(!(i %in% rownames(otu.p))) {
            otus.ex<-unique(c(otus.ex, i))
        }
    } 
    if (length(otus.ex) != 0) {
        print("The following selected otus in your list were not in the otu table:")
        otus.ex    
    }


    # combine otu with taxonomy
    rownames(otu.p.sel) <- paste(rownames(otu.p.sel), otu.p.sel[[rank_name]], sep=":")
    rownames(otu.sel) <- paste(rownames(otu.sel), otu.sel[[rank_name]], sep=":")

    # transpose OTU table
    otu.p.sel.t<-as.data.frame(t(otu.p.sel[, -ncol(otu.p.sel)]))
    otu.p.sel.tax<-otu.p.sel[[rank_name]]
    
    otu.sel.t<-as.data.frame(t(otu.sel[, -ncol(otu.sel)]))
    otu.sel.tax<-otu.sel[[rank_name]]
    
    if(identical(rownames(otu.p.sel.t), rownames(meta))) {
        otu.p.sel.fac <- cbind(otu.p.sel.t, meta[[meta.factor]])
        names(otu.p.sel.fac)[ncol(otu.p.sel.fac)]<-meta.factor
        otu.sel.fac <- cbind(otu.sel.t, meta[[meta.factor]])
        names(otu.sel.fac)[ncol(otu.sel.fac)]<-meta.factor
    } else {
        stop("samples not identical in otu and metadata")
    }
       
    # melt by sample
    #if (!require("reshape2")) {
   # stop("package 'reshape2' is required to use this function")
 # }
    
    otu.sel.fac.m<-melt(cbind(otu.sel.fac, sample=rownames(otu.sel.fac)), 
                             tax=otu.sel.fac$tax, variable.name="OTU", value.name="Count")

    # boxplot or barplot
    if(boxplot) {
        # boxplot
        otu.p.sel.fac.m<-melt(cbind(otu.p.sel.fac, sample=rownames(otu.p.sel.fac)), 
                                     variable.name="OTU", value.name="RA")
        otu.p.sel.fac.m<-cbind(otu.p.sel.fac.m, do.call(rbind,
                            strsplit(x=as.character(otu.p.sel.fac.m$OTU),split=":")))
        names(otu.p.sel.fac.m)[ncol(otu.p.sel.fac.m)-1] <- "otu"
        names(otu.p.sel.fac.m)[ncol(otu.p.sel.fac.m)] <- "tax"
        otu.p.sel.fac.m <- otu.p.sel.fac.m[order(otu.p.sel.fac.m$tax, otu.p.sel.fac.m$otu), ]

    
    } else {
        # barplot
        otu.sel.fac.agg = aggregate(otu.sel.t, by=list(meta[[meta.factor]]), FUN=sum)
        otu.sel.fac.agg <-cbind(otu.sel.fac.agg[,1], decostand(otu.sel.fac.agg[, -1], MARGIN=2, "total"))
        names(otu.sel.fac.agg)[1]<-meta.factor
        # label counts as "RA" to be consistent with previous df
        otu.sel.fac.agg.m <- melt(otu.sel.fac.agg, , is.vars=c(meta.factor), value.name="RA", variable.name="OTU")
        otu.sel.fac.agg.m<-cbind(otu.sel.fac.agg.m, do.call(rbind,strsplit(x=as.character(otu.sel.fac.agg.m$OTU),split=":")))
        names(otu.sel.fac.agg.m)[ncol(otu.sel.fac.agg.m)-1] <- "otu"
        names(otu.sel.fac.agg.m)[ncol(otu.sel.fac.agg.m)] <- "tax"
        otu.sel.fac.agg.m <- otu.sel.fac.agg.m[order(otu.sel.fac.agg.m$tax, otu.sel.fac.agg.m$otu), ]
        otu.p.sel.fac.m <- otu.sel.fac.agg.m
    }
        
    #calculate total count of each Taxon
    total= sapply(levels(otu.sel.fac.m$OTU), function(x){sum(otu.sel.fac.m$Count[otu.sel.fac.m$OTU==x])}, USE.NAMES=F)
    total.count<-data.frame(OTU=levels(otu.sel.fac.m$OTU),total)
    # split otu&tax
    total.count<-cbind(total.count, do.call(rbind,strsplit(x=as.character(total.count$OTU),split=":")))
    names(total.count)[ncol(total.count)-1] <- "otu"
    names(total.count)[ncol(total.count)] <- "tax"
    #total.count <- total.count[order(total.count$tax, total.count$OTU, total.count$total, decreasing=TRUE), ]
        
    #return(list(otu.sel.fac.m, otu.p.sel.fac.m, total.count))

    # reorder levels based on order of appearance in total.count, 
      otu.p.sel.fac.m$OTU <- factor(as.character(otu.p.sel.fac.m$OTU), levels=total.count$OTU, ordered=TRUE)
    

     # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  # the distribution of each taxon in each category of meta.factor
  #  if (!require("scales")) {
  #  stop("package 'scales' is required to use this function")
  #}
  #  if (!require("RColorBrewer")) {
  #  stop("package 'RColorBrewer' is required to use this function")
  #}
  #  if (!require("grid")) {
  #  stop("package 'grid' is required to use this function")
  #}
     
    if (isTRUE(log)) {
        otu.p.sel.fac.m$RA <- log(otu.p.sel.fac.m$RA+1)
        ylab<-"log (Relative Abundance + 1)"
    } else {
        otu.p.sel.fac.m$RA <- otu.p.sel.fac.m$RA
        ylab <- "Relative Abundance"
    }

    if(boxplot) {
    
          p <- ggplot(otu.p.sel.fac.m, aes_string(x="OTU", y="RA",fill=meta.factor)) + 
             geom_boxplot() 
    } else {
        p <- ggplot(otu.p.sel.fac.m, aes_string(x="OTU", y="RA",fill=meta.factor)) + 
           geom_bar(colour="white",stat="identity",position="fill")
    }
    
    p <- p + scale_fill_brewer(name=meta.factor, type="div", palette="Set3") +  
        scale_x_discrete(labels = paste(total.count$otu, "\n", total.count$tax,":",total.count$total, sep="")) + 
        theme(legend.key.size=unit(6,"mm"), axis.text.y=element_text(size=7),legend.position="right") + 
        theme(legend.text = element_text(size=10)) + xlab("OTU") +
        ylab(ylab) + 
        coord_flip()

    if (isTRUE(log)) {
        p <- p
    } else {
        p <- p + scale_y_continuous(labels=percent_format())
    }
    
    print(p)

     
  if (save) { dev.off() }
  
            invisible()
  
}


