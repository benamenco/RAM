percent.classified <- function(data, ranks=c("f","g")) {

  valid.OTU(data)
  
  rk.classified <- list()
  for ( rk in ranks ) {
      .valid.rank(rk)
      rank_name<- .get.rank.name(rk)
      rank_pat<-.get.rank.pat(rk)
      blacklist <- paste(.blacklist(rk), "|", rank_pat, "$", sep="")
  
      # number of rows in get.rank gives us the number of OTUs classified at that
      # rank; divide by total number of OTUs and multiply to get %
      # return(100 * dim(get.rank(otu1=data, rank=rank))[1] / dim(data)[1])
    
      data_rank<-get.rank(otu1=data, rank=rk)
          
      if(length(grep(blacklist, data_rank$taxonomy)) == 0) {
          data_ided <-data_rank
      } else {      
          data_black<-data_rank[grepl(blacklist, data_rank$taxonomy),]
          data_ided<-data_rank[-which(rownames(data_rank) %in% rownames(data_black)),]
      }
      name <- paste("%_OTU_classified_at_", rank_name, "_level:", sep="")
      rk.classified[[name]] <- (100 * nrow(data_ided) / nrow(data))
    }   
  if (length(ranks) == 1) { 
      return(rk.classified[[1]])
  } else {
      return(rk.classified)
  }
}


OTU.recap <- function(otu1, otu2=NULL, ranks=c("g"), labels=NULL) {
  if(is.null(otu2)) {
    num.otus<-1
  } else {
    num.otus<-2
  }
  
  if (!is.null(labels)) {
    if(!(length(labels) ==num.otus)) {
      stop("Error: please provide same number of labels for each otu")
    } else {
      labels=labels
    }
  } else {
    labels=1:2
  }
  
  if (length(ranks)==0) {
    stop("Error: please provide taxonomic ranks to summarize, e.g. ranks=c('p', 'o', 'c', 'f', 'g', 's')")
  }

  names<-c("total_otu", "total_seq", "singleton", "percent_singleton", "percent_singleton_seq", "otu", "seq", "identified_otu", "identified_seq", "percent_identified_otu_of_total", "percent_identified_seq_of_total", "percent_identified_otu_at_rank", "percent_identified_seq_at_rank", "identified_tax")

  len <- length(names)
  len.ranks<-length(ranks)  
  
  df<-list()
  index<-1 # index for first otu
  
  otu_tax<-list()
  for (otu in list(otu1, otu2)){
    # exit if otu2 is null
    if (is.null(otu)) { break }
    
    valid.OTU(otu)
    total_OTUs_num <- nrow(otu)
    total_seqs <-sum(otu[, -ncol(otu)])

    # singletons
    singleton <- which(rowSums(otu[, -ncol(otu)]) == 1)

    if( length(singleton) ==0 ) {
        warning("NO singletons (otus that being observed only once)")
	singleton_OTUs_num <- 0
	singleton_OTUs_percent <- 0
	singleton_seqs_percent <- 0
    } else {

      singleton_OTUs_num <-  nrow(otu[singleton, ])
      singleton_OTUs_percent <- round(100*nrow(otu[rowSums(otu[, -ncol(otu)])==1,])/nrow(otu), digits=3)
      singleton_seqs_percent <-  round(100*nrow(otu[rowSums(otu[, -ncol(otu)])==1,])/total_seqs, digits=3)
    }

      df.names=paste(labels[index], names, sep="_") 

      df[[index]]<-data.frame(matrix(vector(), 0, len, dimnames=list(c(), df.names)), stringsAsFactors=F)    
    
    tax<-list()
    for (i in ranks) {
      .valid.rank(i)
      rank<- .get.rank.name(i)
      rank_pat<-.get.rank.pat(i)
      blacklist <- paste(.blacklist(i), "|", rank_pat, "$", sep="")
      
      # All otus with rank_pat, e.g. f__      
      rank_otu <- otu[grep(rank_pat, otu$taxonomy),]
      if ( nrow(rank_otu) == 0 ) {
        rank_otu_num <- 0
        rank_otu_seq <- 0
      } else {
        rank_otu_num <- nrow(rank_otu)
        rank_otu_seq <- sum(rank_otu[, -ncol(rank_otu)])
      }


      # only identified otus at the rank
      rank_otu_ided <- get.rank(otu1=otu, rank=i)
      col.ided <- ncol(rank_otu_ided)

      # split the taxonomy column
      tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species") 
      rank_otu_ided_splitup <- col.splitup(rank_otu_ided, col="taxonomy", 
                                    sep=";", drop=TRUE, names=tax.classes)
      
      
      if ( (unique(rank_otu_ided_splitup[[rank]])) != "" ) {
          # identified taxa at the rank  
          rank_otu_ided_splitup.rank <- cbind(rank_otu_ided[, 1:(col.ided-1)], 
                                           rank_otu_ided_splitup[[rank]])
          names(rank_otu_ided_splitup.rank)[ncol(rank_otu_ided_splitup.rank)] <- rank
          tax[[i]] = unique(rank_otu_ided_splitup.rank[[rank]])
          rank_tax_ided_num = length(tax[[i]])
     
          rank_otu_ided_num <- nrow(rank_otu_ided)
          rank_otu_ided_seq <- sum(rank_otu_ided[, -ncol(rank_otu_ided)])
          rank_otu_ided_num_percent <- round(100*nrow(rank_otu_ided)/nrow(rank_otu), digits=3)       
          rank_otu_ided_seq_percent <- round(100*sum(rank_otu_ided[, -ncol(rank_otu_ided)]) / sum(rank_otu[, -ncol(rank_otu)]),  digits=3)
          percent_identified_otu_of_total <- round(100*nrow(rank_otu_ided)/nrow(otu), digits=3)
          percent_identified_seq_of_total <- round(100*sum(rank_otu_ided[, -ncol(rank_otu_ided)]) / sum(otu[, -ncol(otu)]),  digits=3)

         
      } else {
	  warning(paste("no otus being identified at ", rank, " level"))
	  rank_otu_ided_num <- 0
          rank_otu_ided_seq <- 0
          rank_otu_ided_num_percent <- 0       
          rank_otu_ided_seq_percent <- 0
          percent_identified_otu_of_total <- 0
          percent_identified_seq_of_total <- 0	
          tax[[i]] = ""
          rank_tax_ided_num = 0
      }

      df[[index]][rank,] <- c(total_OTUs_num, total_seqs, singleton_OTUs_num, singleton_OTUs_percent, singleton_seqs_percent, rank_otu_num, rank_otu_seq, rank_otu_ided_num, rank_otu_ided_seq,  percent_identified_otu_of_total, percent_identified_seq_of_total,  rank_otu_ided_num_percent, rank_otu_ided_seq_percent, rank_tax_ided_num)
    }
    otu_tax[[index]] <- tax
    index<-index+1
  }
  
  if ( num.otus == 1 ) { 
	return(df[[1]])
  } 

  if ( num.otus == 2 ) {
    df_tax_diff<-data.frame(matrix(vector(), len.ranks, 4, dimnames=list(c(), c(paste("num_taxa_only_in_", labels[1], sep=""), paste("num_taxa_only_in_", labels[2], sep=""), paste("taxa_only_in_", labels[1], sep=""), paste("taxa_only_in_", labels[2], sep="")))), stringsAsFactors=F)  
    for (i in ranks) {
      rank<- .get.rank.name(i)
      pos <- which(i==ranks) 
      rownames(df_tax_diff)[pos] <- rank
      
      # find what taxa in otu1 not in otu2
      rank_tax_num_only_in_1<-length(setdiff(otu_tax[[1]][[i]], otu_tax[[2]][[i]]))
      rank_tax_num_only_in_2<-length(setdiff(otu_tax[[2]][[i]], otu_tax[[1]][[i]]))
      tax_only_in_1<-paste(setdiff(otu_tax[[1]][[i]], otu_tax[[2]][[i]]), collapse=",")
      tax_only_in_2<-paste(setdiff(otu_tax[[2]][[i]], otu_tax[[1]][[i]]), collapse=",")

      df_tax_diff[rank, 1] <- rank_tax_num_only_in_1
      df_tax_diff[rank, 2] <- rank_tax_num_only_in_2
      df_tax_diff[rank, 3] <- tax_only_in_1
      df_tax_diff[rank, 4] <- tax_only_in_2
    }
    return(list(df[[1]], df[[2]], df_tax_diff))
  } 

}

OTU.diversity<-function(meta, otu1=otu1, otu2=NULL, labels=NULL) {
  if(is.null(otu2)) {
    num.otus<-1
  } else {
    num.otus<-2
  }
  
  if (!is.null(labels)) {
    if(!(length(labels) ==num.otus)) {
      stop("Error: please provide same number of labels for each otu")
    } else {
      labels=labels
    }
  } else {
    labels=1:2
  }
        
  index<-1

  for (otu in list(otu1, otu2)){
      # exit if otu2 is null
      if (is.null(otu)) { break }
      valid.OTU(otu)
      .valid.meta(otu, otu2=NULL,meta)  
    
      otu.t<-transpose.OTU(otu)
    
    # add columns for different diversity indices
    meta[[paste("spec", labels[index], sep="_")]]<-specnumber(otu.t)
    meta[[paste("sim", labels[index], sep="_")]]<-diversity(otu.t, index="simpson", MARGIN=1)
    meta[[paste("invsim", labels[index], sep="_")]]<-diversity(otu.t, index="invsimpson", MARGIN=1)
    meta[[paste("shan", labels[index], sep="_")]]<-diversity(otu.t, index="shannon", MARGIN=1)
    meta[[paste("sim_even", labels[index], sep="_")]]<-evenness(otu)["otu1",]
    meta[[paste("shan_even", labels[index], sep="_")]]<-evenness(otu, index="shannon")["otu1",]
    meta[[paste("sim_trudiv", labels[index], sep="_")]]<-true.diversity(otu)["otu1",]
    meta[[paste("shan_trudiv", labels[index], sep="_")]]<-true.diversity(otu, index="shannon")["otu1",]
    meta[[paste("chao", labels[index], sep="_")]]<-estimateR(otu.t)["S.chao1",]
    meta[[paste("ACE", labels[index], sep="_")]]<-estimateR(otu.t)["S.ACE",]
    
    index<-index+1
  }
  
  return (meta)
}


