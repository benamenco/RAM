# These functions all deal with creating plots of given data.

######################
# plot.top functions #
######################

### NOTE: since the change in tax.abund 
# (commit 2a420ad061ff9478fc7765d10bd63c2717f614fb), the group.top.[number|percent]
# functions have gotten very messy, and would benefit from a cleanup.

# the current control flow is as follows:
# user calls group.top.[number|percent] -> .top.samples.plot -> 
# top.samples.[ggplot2|base]

# internal function (hidden from user; they call plot.top.percent or plot.top.number)
.top.samples.plot <- function(otu1, otu2=NULL, top=10, drop.unclassified, 
                              labels, file=NULL, ext=NULL, height=8, width=16, 
                              bw=FALSE, ggplot2=TRUE, mode) {
  valid.OTU(otu1, otu2)
  .valid.plot.settings(file, ext)
  
  if (is.null(otu2)) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  # set the function to select the top samples
  # mode, top, and drop.unclassified inherit from this function 
  # (.top.samples.plot) call, data and rank are determined from when top.func 
  # is called later
  top.func <- function(data, rank){tax.abund(data, top=top, rank=rank, 
                                              count=FALSE, mode=mode,
                                              drop.unclassified=drop.unclassified)}
  
  # set the appropriate titles
  if (mode == "number") {
    
    top.title <- paste("Relative Abundance of Top", top, 
                       "Taxon Groups at Five Taxonomic Ranks", sep=" ")
    
  } else if (mode == "percent") {
    top.title <- paste("Relative Abundance of Taxon Groups Above ", top, 
                       "% Relative Abundance at Five Taxonomic Ranks", sep="")
  }
  
  # call the appropriate plotting function
  if (ggplot2) {
    .top.samples.ggplot2(otu1, otu2, top, labels, file, ext, height, width, bw, 
                         top.func, top.title)
  } else {
    .top.samples.base(otu1, otu2, top, labels, file, ext, height, width, bw, 
                      top.func, top.title)
  }
}

.top.samples.ggplot2 <- function(otu1, otu2=NULL, top=10, labels, file=NULL, ext=NULL, 
                                 height=8, width=16, bw=FALSE, top.func, top.title) {
  
  valid.OTU(otu1, otu2)
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  ranks <- c("phylum", "class", "order", "family", "genus")
  
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  # set up our list
  rows <- vector(mode="list", length=length(ranks) * num.otus)
  
  index <- 1
  for (OTU in list(otu1, otu2)){
    # exit if otu2 is null
    if (is.null(OTU)) { break }
    
    # get the number of samples
    num.samples <- dim(OTU)[2] - 1
    
    for (rank in ranks) {
      # choose otu1/otu2 title appropriately
      if (index <= 5) {
        reg <- labels[1]
      } else {
        reg <- labels[2]
      }
      
      # get the data, add region/rank information
      top.otus <- top.func(OTU, rank=rank)
      top.otus <- cbind(top.otus, Region=rep(reg, times=num.samples),
                        Rank=rep(.capitalize(rank), times=num.samples))
      
      # add to our list after melting
      rows[[index]] <- melt(top.otus, id.vars=c("Region", "Rank"), 
                            value.name="RA", variable.name="OTU") 
      index <- index + 1
    }
  }
  
  data <- do.call(rbind, rows)
  
  # get a colour palette with colours for each OTU
  cols.needed <- length(unique(data$OTU))
  if (bw) {
    colours <- rep("white", times=cols.needed)
  } else {
    if (cols.needed > 12) {
      colours <- colorRampPalette(brewer.pal(12, "Set3"))(cols.needed)
    } else {
      colours <- brewer.pal(12, "Set3")
    }
  }
  
  p <- ggplot(data, aes_string(x="OTU", y="RA", fill="OTU")) + 
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=colours) +
    scale_y_continuous(labels = percent_format()) +
    xlab("Taxon") + ylab("Relative Abundance") + 
    ggtitle(top.title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position="none",
          panel.background = element_rect(fill="grey90"),
          panel.grid.major.x = element_line(colour="white"),
          panel.grid.major.y = element_line(colour="white")) +
    facet_grid(Region ~ Rank, scales="free_x")
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}

.top.samples.base <- function(otu1, otu2=NULL, top=10, labels, file=NULL, ext=NULL, 
                              height=8, width=16, bw=FALSE, top.func, top.title) {
  single.otu <- is.null(otu2)
  .valid.plot.settings(file, ext)
  
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  save <- !is.null(file)
  
  if (save) { .get.dev(file, ext, height=height, width=width) }
  
  # backup plot settings and restore afterwards
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  # configure plot settings
  par(oma=c(4, 0, 2, 0), mar=c(8, 2, 1, 1), las=2)
  
  # setup our layout 
  layout(matrix(1:(5 * num.otus), ncol=5, byrow=TRUE))
  
  for (elem in list(otu1, otu2)) {
    # stop if otu2 is null
    if (is.null(elem)) { break }
    
    elem.ranks <- vector(mode="list", length=5)
    
    # call get.rank for each taxonomic classification
    index <- 1
    for (i in c("p", "c", "o", "f", "g")) {
      elem.ranks[[index]] <- top.func(elem, rank=i)
      index <- index + 1
    }
    
    # get color palette
    if (bw) {
      cols <- grey.colors(top)
    } else {
      # the max pallette size is 12; if length > 12 we get the palette with 12
      # colors (which will be recycled by boxplot, which is fine)
      cols <- suppressWarnings(brewer.pal(top, "Set3"))
    }
    
    tax.names <- c("Phylum", "Class", "Order", "Family", "Genus")
    
    # plot each rank
    for (i in 1:length(tax.names)) {
      # plot the data
      boxplot(elem.ranks[[i]], col=cols, notch=FALSE, cex.axis=0.8)
      # add a legend
      if (identical(elem, otu1)) {
        region <- labels[1]
      } else {
        region <- labels[2]
      }
      legend("topright", legend=paste0(region, " - ", tax.names[i]))
    }
  }
  
  # print title
  title(main=top.title, outer=TRUE)
  
  if (save) { dev.off() }
  
  invisible()
}

# group.top.percent and group.top.number are what the user call; the other top.samples.X
# functions are used interally to produce the graph.
group.top.percent <- function(otu1, otu2=NULL, top=10, drop.unclassified=FALSE, 
                              labels=c("ITS1", "ITS2"), 
                              file=NULL, ext=NULL, height=8, width=16, 
                              bw=FALSE, ggplot2=TRUE) {
  valid.OTU(otu1, otu2)
  # quickly get the numer of OTUs, !is.null will be coerced to 0 or 1
  .valid.labels(1 + !is.null(otu2), labels)
  
  .top.samples.plot(otu1, otu2, top=top, labels=labels, 
                    drop.unclassified=drop.unclassified,
                    file=file, ext=ext, height=height, width=width,
                    bw=bw, ggplot2=ggplot2, mode="percent")
}

group.top.number <- function(otu1, otu2=NULL, top=10, drop.unclassified=FALSE,
                             labels=c("ITS1", "ITS2"), 
                             file=NULL, ext=NULL, height=8, width=16, 
                             bw=FALSE, ggplot2=TRUE) {
  valid.OTU(otu1, otu2)
  # quickly get the numer of OTUs, !is.null will be coerced to 0 or 1
  .valid.labels(1 + !is.null(otu2), labels)
  
  .top.samples.plot(otu1, otu2, top=top, labels=labels, 
                    drop.unclassified=drop.unclassified,
                    file=file, ext=ext, height=height, width=width, 
                    bw=bw, ggplot2=ggplot2, mode="number")
}

group.abundance <- function(otu1, otu2=NULL, rank, 
                            top=NULL, count=FALSE, drop.unclassified=FALSE,
                            file=NULL, ext=NULL, labels=c("ITS1", "ITS2"),
                            height=8, width=16, bw=FALSE, ggplot2=TRUE) {
  
  .valid.plot.settings(file, ext)
  .valid.rank(rank)
  
  if (is.null(otu2)) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  .valid.labels(num.otus, labels)
  
  if (ggplot2) {
    .abundance.ggplot2(otu1, otu2, rank, top, count, drop.unclassified, file, ext,
                       labels, height, width, bw)
  } else {
    .abundance.base(otu1, otu2, rank, top, count, drop.unclassified, file, ext,
                    labels, height, width, bw)
  }
}

.abundance.base <- function(otu1, otu2=NULL, rank, 
                            top=NULL, count=FALSE, drop.unclassified=FALSE,
                            file=NULL, ext=NULL, labels, 
                            height=8, width=9, bw=FALSE) {
  valid.OTU(otu1, otu2)
  .valid.rank(rank)
  .valid.plot.settings(file, ext)
  
  save <- !is.null(file)
  
  single.otu <- is.null(otu2)
  
  if (single.otu) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  # set up the appropriate device if saving
  if (save) {
    .get.dev(file, ext, height, width)
  }
  
  # backup plot settings and restore afterwards
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  # set up the layout: one shared pane up top (for the title),
  # one pane for each graph and one for its legend (if applicable)
  #mat <- matrix(c(rep(1, times=num.otus), 2:((2 * num.otus) + 1)), ncol=num.otus,
  #              byrow=TRUE)
  
  mat <- matrix(c(1,1,2,4,3,5), ncol=2, byrow=TRUE)
  
  if (num.otus == 1) {
    lmat <- matrix(1:3, ncol=1)
  } else {
    lmat <- matrix(c(1:3, c(1, 4, 5)), ncol=2)
  }
  
  layout(lmat, heights=c(1, 8, 2))
  
  pretty.rank <- .get.rank(.get.rank.ind(rank), pretty=TRUE)
  
  if (count) {
    top.title <- paste("Counts of Taxonomic Groups at", pretty.rank, "Level") 
  } else {
    top.title <- paste("Relative Abundance of Taxonomic Groups at", pretty.rank,
                       "Level")
  }
  
  par(mar=c(0, 0, 4, 0))
  plot.new()
  title(main=top.title, cex.main=1.5)
  
  index <- 1
  
  for (elem in list(otu1, otu2)) {
    # stop if not given otu2
    if (is.null(elem)) { break }
    
    abund <- tax.abund(elem, rank=rank, drop.unclassified=drop.unclassified,
                         top=top, count=count)
    
    sample.totals <- c(rowSums(abund), 0)
    breaks <- pretty(sample.totals)
    
    num.taxa <- dim(abund)[2]
    
    # set the appropriate label
    if (count) {
      y.label <- "Count"
    } else {
      y.label <- "Relative Abundance"
    }
    
    if (bw) {
      cols <- grey.colors(n=num.taxa, end=0.7)
    } else {
      # generate palette of colours, warn the user if too many OTUs
      cols <- suppressWarnings(brewer.pal(num.taxa, "Set3"))
      
      if (num.taxa > 12) {
        warning("this colour palette only has 12 distinct colours, and more than 12 taxon groups have been provided. Some colours are being recycled in your graph (be careful!)")
      }
    }
    
    # the barplot function is expecting the transpose
    abund <- t(abund)
    
    # we add the axes manually with our breaks variable later
    barplot.args <- list(abund, beside=FALSE, names.arg=colnames(abund),
                         cex.names=0.7, col=cols, las=2, axes=FALSE,
                         ylim=c(0, max(breaks)))
    
    leg.args <- list(x="center", legend=rownames(abund), cex=0.7, horiz=FALSE, 
                     ncol=3, fill=cols)
    
    # add shading to the plot
    if (bw) {
      dens <- rep(c(18, 9), length.out=num.taxa)
      ang <- c(30, 60, 90, 120, 150)
      
      barplot.args <- c(barplot.args, density=list(dens), angle=list(ang))
      leg.args <- c(leg.args, density=list(3 * dens), angle=list(ang))
    }
    
    par(mar = c(5.1, 5.1, 4.1, 2.1))
    do.call(barplot, barplot.args)
    
    title(ylab=y.label, line=3.75)
    title(xlab="Samples", line=4.2)
    title(main=labels[index])
    
    
    axis(side=2, at=breaks, las=2)
    
    
    # set margins for legend
    par(mar=c(0, 0, 0, 0))
    plot.new()
    do.call(legend, leg.args)
    index <- index + 1
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}

.abundance.ggplot2 <- function(otu1, otu2=NULL, rank, 
                               top=NULL, count=FALSE, drop.unclassified=FALSE,
                               file=NULL, ext=NULL, labels, 
                               height=8, width=16, bw=FALSE) {
  # validate inputs
  valid.OTU(otu1, otu2)
  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  single.otu <- is.null(otu2)
  
  taxa <- vector(mode="list", length=2)
  index <- 1
  
  for (elem in list(otu1, otu2)) {
    if (is.null(elem)) {
      break
    }
    # get the groups
    elem.tax <- tax.abund(elem, rank=rank, drop.unclassified=drop.unclassified,
                          count=count, top=top)
    # melt by Sample
    elem.tax <- melt(cbind(elem.tax, Sample=rownames(elem.tax)), id.vars="Sample",
                     value.name="RA", variable.name="Rank")
    # bind together with appropriate label
    elem.tax <- cbind(elem.tax, Region=labels[index])
    
    taxa[[index]] <- elem.tax
    index <- index + 1
  }
  
  all.taxa <- rbind(taxa[[1]], taxa[[2]])
  
  # reorder levels based on order of appearance in table
  all.taxa$Sample <- ordered(all.taxa$Sample, levels=unique(all.taxa$Sample))
  
  if (count) {
    title <- paste("Counts of Taxonomic Groups at", 
                   .get.rank(.get.rank.ind(rank), pretty=TRUE),
                   "Level")
  } else {
    title <- paste("Relative Abundance of Taxonomic Groups at", 
                   .get.rank(.get.rank.ind(rank), pretty=TRUE),
                   "Level")
  }
  
  # we need to use aes_string to pass CRAN check; see 
  # http://goo.gl/JxgZ9u
  p <- ggplot(all.taxa, aes_string(x="Sample", y="RA", fill="Rank")) + 
    geom_bar(position="stack", stat="identity") +
    #coord_flip() +
    xlab("Samples") +
    ggtitle(title) +
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major.x = element_blank())
  
  if (!count) {
    p <- p +  scale_y_continuous(labels = percent_format()) +
      ylab("Relative Abundance")
  } else {
    p <- p + ylab("Count")
  }
  
  if (!single.otu) { # for multiple OTUs, wrap based on region
    p <- p + facet_wrap(~Region, scales="free_x")
  }
  
  if (bw) { # for black/white plots
    warning("the ggplot2 package used to create this graph cannot handle patterning for bar plots; greyscale shading is being used.")
    p <- p + scale_fill_grey(guide = guide_legend(direction="horizontal", ncol=5)) +
      theme_bw() + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
  } else { # use actual colours
    
    cols.needed <- length(unique(all.taxa$Rank))
    
    if (cols.needed > 12) {
      # palette max is 12; so if we have more than 6 entries for otu1/2, we 
      # need to manually construct a palette to use
      
      col.func <- colorRampPalette(brewer.pal(12, "Set3"))
      p <- p + scale_fill_manual(values=col.func(cols.needed), 
                                 guide=guide_legend(direction="horizontal", ncol=5))
    } else {
      # otherwise, we can just use the Set3 palette from RColorBrewer
      p <- p + scale_fill_brewer(palette="Set3",
                                 guide=guide_legend(direction="horizontal", ncol=5))
    }
  }
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}

pcoa.plot <- function(data, meta, factors, rank, stand.method=NULL,
                      dist.method="morisita", sample.labels=TRUE,
                      top=20, ellipse=FALSE, main=NULL, file=NULL, ext=NULL,
                      height=8, width=10, ggplot2=TRUE, bw=FALSE) {
  
  valid.OTU(data)
  .valid.rank(rank)
  .valid.meta(otu1=data, meta=meta)
  num.facs <- length(factors)
  
  if (!is.numeric(top) || length(top) != 1L || top < 0) {
    stop("'top' must be a numeric vector of length 1 and have a value > 0.")
  }
  
  if (!is.logical(sample.labels) || length(sample.labels) != 1L || 
        is.na(sample.labels)) {
    stop("'sample.labels' must be a logical vector of length 1.")
  }
  
  if (!is.logical(bw) || length(bw) != 1L || is.na(bw)) {
    stop("'bw' must be a logical vector of length 1.")
  }
  
  # we have to check with identical, because in R, TRUE == 1 but 
  # identical(TRUE, 1) is FALSE
  if (!any(identical(ellipse, 1), identical(ellipse, 2), 
           identical(ellipse, FALSE))) {
    stop("'ellipse' should be either 1, 2, or FALSE; see ?pcoa.plot for details.")
  }
  
  if (ellipse > num.facs) {
    stop("argument 'ellipse' cannot be greater than the number of factors.")
  }
  
  if (ellipse && ggplot2) {
    warning("drawing the ellipses for groups is not currently supported when 'ggplot2=TRUE'.")
  }
  
  # get the metadata factors
  meta.factors <- .valid.factors(meta, factors, min.factors=1, max.factors=2)
  
  # if we're plotting 2 factors in b&w, we need to 'cross' the two factors and 
  # assign a unique symbol to each
  if (bw && num.facs >= 2) {
    
    # ggplot complains if we include /, -, " ", etc. in names
    if (ggplot2) {separator <- "_"} else {separator <- "/"}
    
    cross.name <- paste(names(meta.factors), collapse=separator)
    meta.factors <- cbind(paste(meta.factors[ ,1], meta.factors[ ,2], 
                                sep=separator),
                          meta.factors)
    
    names(meta.factors)[1] <- cross.name
  }
  
  for (column in colnames(meta.factors)) {
    if (!is.factor(meta.factors[ ,column])) {
      warning(paste("column", column,
                    "from 'factors' is not a factor; attempting to coerce now", 
                    "(see ?RAM.factors for help)."))
      
      meta.factors[ ,column] <- as.factor(meta.factors[ ,column])
    }
    
    if (length(levels(meta.factors[ ,column])) > 9) {
      warning(paste("there are more than 9 levels in column", column,
              "from 'factor'; only 9 will be shown."))
    }
    
    # to plot ellipses, we need more than 2 counts for each level
    if (any(summary(meta.factors[ ,column]) <= 2) && ellipse) {
      warning(paste("column", column, "from 'factor' has less than two",
                    "observations for a level, this prevents ellipses from",
                    "being plotted."))
    }
  }
  
  otu.t <- transpose.OTU(data)
  abund <- tax.abund(data, rank=rank, drop.unclassified=TRUE)
  
  if (!is.null(stand.method)) {
    abund <- decostand(abund, stand.method)
  }
  
  dists <- vegdist(abund, method=dist.method)
  
  k.max <- dim(otu.t)[1] - 1
  
 # if (!require("labdsv")) {
 #     stop("package 'labdsv' is required to use this function")
 # }
  
  pcoa <- suppressWarnings(pco(dists, k.max))
  
  # if we don't get at least two axes of ordination, throw error
  if (dim(pcoa$points)[2] < 2) {
    stop("less than two axes of ordination were calculated. Please try again with different values for 'stand.method' and/or 'dist.method'.")
  }
  
  sp.scores <- wascores(pcoa$points, abund)
  
  if (ggplot2) {
    .pcoa.ggplot2(abund, pcoa, rank, sp.scores, meta.factors, sample.labels, 
                  top, ellipse, main, file, ext, height, width, bw)
  } else {
    .pcoa.base(abund, pcoa, rank, sp.scores, meta.factors, sample.labels, top, 
               ellipse, main, file, ext, height, width, bw)
  }
}

.pcoa.ggplot2 <- function(abund, pcoa, rank, sp.scores, meta.factors, 
                          sample.labels, top, ellipse, main, file, ext, height, width,
                          bw) {
  
  num.facs <- length(meta.factors)
  save <- !is.null(file)
  
  # set up the data frame with pcoa data
  samples <- as.data.frame(cbind(Sample=rownames(pcoa$points), 
                                 X=pcoa$points[ , 1], Y=pcoa$points[ , 2]))
  
  for (i in 1:num.facs) {
    samples <- cbind(samples, meta.factors[[i]])
    names(samples)[i + 3] <- names(meta.factors)[i]
  }

  samples$X <- as.numeric(as.character(samples$X))
  samples$Y <- as.numeric(as.character(samples$Y))
  samples$Sample <- as.character(samples$Sample)
  
  # set up the data frame with taxonomic data
  otus <- as.data.frame(cbind(OTU=rownames(sp.scores), X=sp.scores[ ,1], 
                              Y=sp.scores[ ,2]))
  otus$X <- as.numeric(as.character(otus$X))
  otus$Y <- as.numeric(as.character(otus$Y))
  
  # filter out the top samples
  if (dim(otus)[1] < top) {
    warning(paste("there are less than", top, 
                  "taxon groups at the given rank; plotting them all."))
    
    top <- dim(sp.scores)[1]
  }
  
  otus <- otus[1:top, ]
  
  # determine aes based on number of meta factors and bw setting
  # (recall the special case when bw=T and num.facs >=2 since we 'crossed' the
  # factors)
  if (num.facs >= 2) {
    if (bw) {
      samples.aes <- aes_string(x="X", y="Y", label="Sample",
                                shape=names(samples)[4])
    } else {
      samples.aes <- aes_string(x="X", y="Y", label="Sample", 
                                shape=names(samples)[4],
                                colour=names(samples)[5])
    }
  } else if (num.facs == 1) {
    samples.aes <- aes_string(x="X", y="Y", label="Sample", shape=names(samples)[4])
  }
  
  num.samples <- length(unique(samples$Sample))
  
  # create a vector for the vertical position of labels; make half the labels 
  # justified up/down on y-axis, the others constant on y-axis
  v.jitter <- sample(c(2.8, 0.5), size=num.samples, replace=TRUE)
  v.jitter <- sapply(v.jitter, 
                     FUN=function(x){if (x != 0.5 && runif(1) < 0.5) {x * -1} else {x}})
  
  # for the labels with no y-axis jitter, add x-axis jitter
  h.jitter <- sapply(v.jitter, 
                     FUN=function(x){if (x != 0.5) {0.5} else {sample(c(-0.4, 1.4), size=1)}})
  
  
  # add points for samples (with labels)
  p <- ggplot(samples, samples.aes) + geom_point(alpha=0.65, size=7)
  
  if (sample.labels) {
    # jitter the sample labels (using vectors from above)
    p <- p + geom_text(size=2, colour="black", vjust=v.jitter, hjust=h.jitter)
  }
    
  if (top != 0) {
    # add taxon groups
    p <- p +
      geom_text(aes_string(x="X", y="Y", label="OTU", colour=NULL, shape=NULL),
                data=otus, size=3, colour="darkgrey", alpha=0.8)
  }
  
  x.lab <- paste0("Axis I (", round(100 * pcoa$eig[1] / sum(pcoa$eig), digits=2), "%)")
  y.lab <- paste0("Axis II (", round(100 * pcoa$eig[2] / sum(pcoa$eig), digits=2), "%)")
  
  #main.title <- paste("PCoA for Top", top, "Taxon Groups at",
  #                    .get.rank(.get.rank.ind(rank), pretty=TRUE),
  #                     "Level")
  if (is.null(main)) {
        main.title = ""
  } else {
      main.title = main
  }
  
  # add titles
  p <- p + ggtitle(main.title) + xlab(x.lab) + ylab(y.lab)
  
  # add theme
  p <- p + theme(panel.background = element_rect(fill="white"),
                 panel.border = element_rect(colour="black", fill="transparent"))

  # add colour
  p <- p + scale_color_brewer(palette="Set1")
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  }
  
  p
}

.pcoa.base <- function(abund, pcoa, rank, sp.scores, meta.factors, 
                       sample.labels, top, ellipse, main, file, ext, height, width,
                       bw) {
  
  .valid.plot.settings(file, ext)
  save <- !is.null(file)
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  num.facs <- min(length(meta.factors), 2)
  
  # ugly hack, if bw=T and we have >=2 factor, we've "crossed" two other factors
  # so manually set num.facs <- 1
  
  if (bw && num.facs >=2) { num.facs <- 1 }
  
  if (save) {
    .get.dev(file, ext, height=height, width=width)
  }
  
  lmat <- matrix(c(rep(1, times=num.facs), 2:(num.facs + 1)), ncol=2)
  layout(lmat, widths=c(9, 1))
  
  # get label names
  x.lab <- paste0("Axis I (", round(100 * pcoa$eig[1] / sum(pcoa$eig), digits=2), "%)")
  y.lab <- paste0("Axis II (", round(100 * pcoa$eig[2] / sum(pcoa$eig), digits=2), "%)")
  # main.title <- paste("PCoA for Top", top, "Taxon Groups at",
  #                    .get.rank(.get.rank.ind(rank), pretty=TRUE),
  #                    "Level")
  
  if (is.null(main)) {
      main.title = ""
  } else {
      main.title = main
  }
  
  plot.args <- list(pcoa$points[ , 1:2], xlab=x.lab, ylab=y.lab, main=main.title,
                    lwd=4, cex=2.5)
  
  par(mar=c(5.1, 4.1, 4.1, 0))
  
  # set up colours for meta factors 
  meta.cols <- vector(length=num.facs, mode="list")
  palettes <- list(brewer.pal(9, "Pastel1"), brewer.pal(9, "Set1"))
  
  for (i in 1:num.facs) {
    
    meta.cols[[i]] <- palettes[[i]][as.numeric(meta.factors[[i]])]
    names(meta.cols[[i]]) <- as.character(meta.factors[[i]])
  }
  
  # the shapes to be used for bw plotting (we pass these numbers to pch)
  shapes <- c(0:2, 5:6, 15:18)
  
  if (num.facs == 2) {
    
      plot.args <- c(plot.args, pch=21, bg=list(meta.cols[[1]]), 
                     col=list(meta.cols[[2]]))
    
  } else { # only one metadata factor
    if (bw) {
      plot.args <- c(plot.args, pch=list(shapes[as.numeric(meta.factors[[1]])]))
      
    } else {
      plot.args <- c(plot.args, pch=16, col=list(meta.cols[[1]]))
    }
  }
  
  do.call(plot, plot.args)
  
  if (ellipse) {
    
    # if bw=TRUE & num.facs >=2, the first column is the "crossed" factor, which
    # we correct for by adding 1
    if (bw && num.facs >= 2) { ellipse <- ellipse + 1}
    
    # add ellipses with group names as centroids
    for (level in unique(meta.factors[[ellipse]])) {
      
      if (bw) {
        ordiellipse(pcoa, meta.factors[[ellipse]], show.groups = level,
                    col="black", label=TRUE)
        
      } else { # plot in colour
        ordiellipse(pcoa, meta.factors[[ellipse]], show.groups = level,
                    # this call gets the appropriate colour palette, then extracts
                    # the colour corresponding to the value of the factor level
                    col=meta.cols[[ellipse]][names(meta.cols[[ellipse]]) == level],
                    label=TRUE)
      }
    }
  }
  
  if (sample.labels) {
  # plot the sample labels
  text(pcoa$points[ , 1:2], labels=rownames(abund), 
       # randomly place the labels; call the function until this looks nice
       pos=sample(4, size=length(rownames(abund)), replace=TRUE), 
       cex=0.6, offset=1)
  }
  
  # filter out the top samples
  if (dim(sp.scores)[1] < top) {
    warning(paste("there are less than", top, 
                  "taxon groups at the given rank; plotting them all."))
    
    top <- dim(sp.scores)[1]
  }
  
  if (top != 0) {
    # plot the taxonomic information
    text(sp.scores[1:top, 1:2, drop=FALSE], labels=rownames(sp.scores)[1:top], 
         col="darkgrey", cex=0.8)
  }
  
  for (i in 1:num.facs) {
    
    # align the legend boxes with the edges of the plot
    if (i == 1) {
      par(mar=c(0, 0, 4.1, 0))
      plot.new()
      dir <- "topleft"
      # pch value for filled dot
      leg.shape <- 16
      
    } else if (i == 2) {
      par(mar=c(5.1, 0, 0, 0))
      plot.new()
      dir <- "bottomleft"
      # pch value for filled dot
      leg.shape <- 21
    }
    
    leg.args <- list(dir, legend=levels(meta.factors[[i]]), 
                     title=names(meta.factors)[i], xpd=NA)
    
    # if black and white, use shapes for legend, otherwise use colour
    if (bw) {
      leg.args <- c(leg.args, pch=list(shapes[1:length(levels(meta.factors[[i]]))]))
      
    } else {
      leg.args <- c(leg.args, pch=leg.shape, bg="white", pt.lwd=3, pt.cex=1.5,
                    col=list(unique(meta.cols[[i]])))
    }
    
    do.call(legend, leg.args)
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}

group.temporal <- function(data, meta, date.col, factors, rank, group, 
                            file=NULL, ext=NULL, height=8, width=12) {
  
  valid.OTU(data)
  .valid.rank(rank)
  .valid.meta(otu1=data, meta=meta)
  save <- !is.null(file)
  
  # check all sample names (ignore the last column; that's taxonomy)
  if (!all(colnames(data)[-dim(data)[2]] %in% rownames(meta))) {
    stop("sample names do not match for OTU table and metadata.")
  }
  
  if (!date.col %in% colnames(meta)) {
    stop("'date.col' was not found in 'meta'.")
  }
  
  if (class(meta[ ,date.col]) != "Date") {
    warning("the date column in metadata must be of type Date; coercing it now. See ?RAM.dates.")
    meta[ ,date.col] <- as.Date(meta[ ,date.col])
  }
  
  meta.factors <- .valid.factors(meta, factors)
  
  # regular apply will coerce the data frame to a matrix, so all columns become 
  # character vectors (infuriating!); so instead lapply and unlist after (since
  # data frames are lists in disguise)
  not.numeric <- unlist(lapply(meta.factors,
                       FUN=function(column) {
                         !is.numeric(column)}))
  
  if (any(not.numeric)) {
    stop((paste("this function only supports numeric variables, the following variables are not numeric:\n",
                paste(names(meta.factors)[not.numeric], collapse=" "))))
  }
  
  dates <- meta[ ,date.col]
  
  meta.factors <- aggregate(meta.factors, by=list(dates), FUN=mean)
  meta.factors <- rename(meta.factors, c("Group.1"="Date"))
  meta.factors <- melt(meta.factors, id.vars="Date", variable.name="Measure",
                       value.name="Value")
  
  xlims <- c(min(meta.factors$Date), max(meta.factors$Date))
  
  meta.plot <- ggplot(meta.factors, aes_string(x="Date", y="Value")) + geom_line() +
               facet_wrap(~Measure, ncol=1, scales="free_y") + 
               ylab("Average Value") + 
               theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     plot.margin=unit(c(1, 1, 0, 0.5), "lines")) +
               xlim(xlims)
     
  
  groups <- .get.tax.group(data, rep(rank, times=length(group)), group)
  abund <- tax.abund(groups, rank=rank)
  
  ordering <- match(rownames(abund), rownames(meta))
  meta.sort <- meta[ordering, ]

  # tax.abund reorders the samples, so be careful when aggregating
  abund.agg <- aggregate(abund, by=list(dates[ordering]), FUN=sum)
  abund.agg <- rename(abund.agg, c("Group.1"="Date"))
  
  abund.agg <- melt(abund.agg, id.vars="Date", variable.name="Group", 
                     value.name="Count")
  
  main.plot <- ggplot(abund.agg, aes_string(x="Date", y="Count", fill="Group")) + 
               geom_bar(stat="identity", position="stack") +
               theme(legend.position="bottom",
                     plot.margin=unit(c(0, 1, 1, 0.5), "lines")) +
               xlim(xlims) + 
               scale_fill_brewer(palette = "Set3")
  
  # much credit due to baptiste (on SO)
  # taken from http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
  gA <- ggplot_gtable(ggplot_build(meta.plot))
  gB <- ggplot_gtable(ggplot_build(main.plot))

  maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
  gA$widths[2:3] <- as.list(maxWidth)
  gB$widths[2:3] <- as.list(maxWidth)
  
  if (save) {
    file <- .ensure.filepath(file, ext)
    .get.dev(file, ext, height=height, width=width)
  }
  
  grid.arrange(gA, gB, ncol=1, heights=c(0.2, 0.8))
  
  if (save) { dev.off() }
  invisible()
}

group.spatial <- function(data, meta, date.col, province.col, rank, group, 
                          breaks="year", file=NULL, ext=NULL, height=8, width=10) {
  valid.OTU(data)
  .valid.rank(rank)
  .valid.meta(otu1=data, meta=meta)
  save <- !is.null(file)
  
  if (!date.col %in% colnames(meta)) {
    stop("'date.col' was not found in 'meta'.")
  }
  
  if (!province.col %in% colnames(meta)) {
    stop("province.col' was not found in 'meta'.")
  }
  
  if (!is.character(meta[ ,province.col])) {
    warning("the 'Province' column of meta is not a character vector; coercing it to one.")
    meta$Province <- as.character(meta[ ,province.col])
  }
  
  meta$Province <- .get.province.id(meta[ ,province.col])
  
  group.len <- length(group)
  groups <- .get.tax.group(data, rep(rank, times=group.len), group)
  abund <- tax.abund(groups, rank=rank)
  
  ordering <- match(rownames(abund), rownames(meta))
  meta.sort <- meta[ordering, ]
  
  meta.sort[ ,date.col] <- as.Date(meta.sort[ ,date.col])
  # note that setting ordered_result = TRUE in cut does NOT order the 
  # levels by actual date value, so we do that manually
  date.factor <- cut(meta.sort[ ,date.col], breaks=breaks)
  date.factor <- ordered(date.factor, sort(levels(date.factor)))
  
  abund.loc <- aggregate(abund, by=list(meta.sort$Province, date.factor), 
                         FUN=sum)
  
  names(abund.loc)[1:2] <- c("id", "Date")
  
  # this is only here to silence a CRAN check note; otherwise it complains 
  # that we reference map.fortify without defining it (even though loading
  # the file below will create an object map.fortify)
  map.fortify.list <- NULL
  # load file map.fortify, which contains a Canadian province map already 
  # fortified (for plotting with ggplot2)
  load(system.file("extdata", "map.fortify.list.RData", package="RAM"))
  
  # we get all missing provinces (which have count 0), and add them to 
  # abund.loc so they show up in our map later (for each time segment)
  missing <- vector(length=length(levels(abund.loc$Date)), mode="list")
  
  dates <- levels(abund.loc$Date)
  for (i in 1:length(dates)) {
    # get all provinces without counts
    missing.provinces <- setdiff(unique(map.fortify.list$map$id),
                                 unique(abund.loc[abund.loc$Date == dates[i], ]$id))
    
    # add a NA count for each group
    group.count <- vector(length=group.len, mode="list")
    group.count[1:group.len] <- NA
    
    # create and rename data frame
    missing[[i]] <- data.frame(missing.provinces, dates[i], group.count)
    # all names after the first two are groups
    names(missing[[i]]) <- names(abund.loc)
  }
  
  abund.loc <- rbind(abund.loc, do.call(rbind, missing))
  abund.loc.melt <- melt(abund.loc, id.vars=c("id", "Date"), variable.name="Group",
                         value.name="Count")
  
  loc.map <- base::merge(map.fortify.list$map, abund.loc.melt, by="id")
  
  # set the fill to the given taxon group
  tax.group <- names(abund.loc)[3]

  p <- ggplot(loc.map, aes_string(x="long", y="lat", map_id="id", fill="Count")) +
       scale_fill_gradient(na.value = "white") +
       geom_map(map=loc.map, colour="black") +
       facet_grid(Group ~ Date) +
       theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
             axis.text.y = element_blank(), axis.title.y = element_blank(),
             axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
             panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
}

# this function simply maps the province codes to the names used in the GeoBase map
.get.province.id <- function(province) {
  
  valid.abbrev <- c("AB", "BC", "MB", "NB", "NL", "NT", "NS", "NU", "ON", "PE",
                    "QC", "SK", "YT")
  
  provinces <- c("ALBERTA", "BRITISH COLUMBIA", "MANITOBA", "NEW BRUNSWICK", 
                 "NEWFOUNDLAND AND LABRADOR", "NORTHWEST TERRITORIES", "NOVA SCOTIA",
                 "NUNAVUT", "ONTARIO", "PRINCE EDWARD ISLAND", "QUEBEC", 
                 "SASKATCHEWAN", "YUKON" )
  
  if (!is.character(province)) {
    warning("'province' should be a character vector; coercing it to one now.")
    province <- as.character(province)
  }
  
  if (!all(province %in% valid.abbrev)) {
    stop("invalid province code. See ?location.formatting for details.")
  }

  # get the name at the correct index
  provinces[match(province, valid.abbrev)]
}

group.indicators <- function(otu1, otu2=NULL, meta, factor, rank,
                            thresholds = c(A=0.85, B=0.8, stat=0.8, p.value=0.05),
                            labels=c("ITS1", "ITS2"), file=NULL, ext=NULL,
                            height=12, width=12) {
  
  # I have seen this discussion: http://yihui.name/en/2014/07/library-vs-require/
  # but I think returning an explanatory error message is worthwhile
  
  if (require("indicspecies")) {
      indicspecies::multipatt
  } else {
    stop("package 'indicspecies' is required to use this function: try 'install.packages('indicspecies')'.")
  }
  
  # we create some temporary text files in this function, so clean up:
  # (this removes all files containing mp_summary in their name in the R temp dir)
  on.exit(unlink(paste0(tempdir(), "/mp_summary*")))
  
  valid.OTU(otu1, otu2)
  .valid.rank(rank)
  .valid.meta(otu1, otu2, meta)
  save <- !is.null(file)
  
  if (is.null(otu2)) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  .valid.labels(num.otus, labels)
  
  meta.factor <- .valid.factors(meta, factor, min.factors = 1, max.factors = 1)

  if (!is.numeric(thresholds) || length(thresholds) != 4L) {
    stop("thresholds must be a numeric vector of length four (see ?indicators.plot for details).")
  }
  
  meta.name <- names(meta.factor)
   
  # to store the data we're about to generate
  rows <- vector(length=num.otus, mode="list")
  
  index <- 1
  for (otu in list(otu1, otu2)) {
    # stop if not given otu2
    if (is.null(otu)) { break }
    
    abund <- tax.abund(otu, rank=rank)
    abund.stand <- decostand(abund, method="total")
    
    mp <- indicspecies::multipatt(abund.stand, meta.factor[[1]], control=how(nperm=999))
     
    # write the summary of the multi-level analysis to a temp file, which we read
    tempfilename <- tempfile("mp_summary", fileext=".txt")
    sink(file=tempfilename)
    summary(mp, indvalcomp=TRUE)
    sink()
    
    indicators <- readLines(tempfilename)
    
    # the summary contains some human-friendly non-data lines we need to strip
    # this regex keeps only lines that have {number number number number},
    # since there are four numerical columns in the data (and none elsewhere)
    # this selects only the numerical data we want
    matches <- grepl("([[:digit:]\\.]+[[:space:]]+){4}", indicators)
    
    if (!any(matches)) {
      stop("no taxon groups were significant with the given parameters; try again with a lower taxon group and/or a different meta.factor.")
    }
    
    indicators <- indicators[matches]
    
    # split the rows at all whitespace, convert to dataframe
    indicators.df <- data.frame(Group=character(), A=numeric(), B=numeric(),
                                stat=numeric(), p.value=numeric(),
                                stringsAsFactors=FALSE)
    
    pieces <- strsplit(indicators, "[[:space:]]+")
    
    for (i in 1:length(pieces)) {
      # the [-6] removes the trailing significance code (which is just a human-friendly
      # display of the p-values)
      indicators.df[i, ] <- pieces[[i]][-6]
    }
    
    for (i in 2:5) {
      indicators.df[ , i] <- as.numeric(indicators.df[ , i])
    }
   
    # get names of all taxon groups above the given thresholds
    keep <- indicators.df[with(indicators.df, 
                               A >= thresholds[1] & B >= thresholds[2] & 
                               stat >= thresholds[3] & p.value <= thresholds[4])
                          , "Group", drop=TRUE]
    
    if (length(keep) == 0L) {
      stop("no taxon groups met all thresholds. Either relax your thresholds, or try again with each otu individually.")
    }
    
    # set up our data frame with all the data ggplot needs
    abund.filtered <- data.frame(Sample=rownames(abund.stand),
                                 meta=meta.factor,
                                 Region=rep(labels[index], times=dim(abund.stand)[1]), 
                                 abund.stand[ , keep, drop=FALSE])
    
    # melt it & store the result
    rows[[index]] <- melt(abund.filtered, id.vars=c("Sample", meta.name, "Region"),
                                                    variable.name="Indicator",
                                                    value.name="Value")
    
    index <- index + 1
  } # end otu for loop
  
  data <- do.call(rbind, rows)
  
  # we use as.formula/paste to create the faceting formula (otherwise we would 
  # have to type Region ~ meta.name, and facet_grid won't evaluate meta.name)
  p <- ggplot(data, aes_string(x="Sample", y="Value", fill="Indicator")) +
       geom_bar(stat="identity") +
       facet_grid(as.formula(paste("Region", "~", meta.name)),
                  scales="free_x", space="free_x") +
       scale_fill_brewer(palette="Set3") +
       ylab("Relative Abundance") +
       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="bottom")
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p
  }
}

sample.locations <- function(otu1, otu2=NULL, meta, factor=NULL, zoom=5,
                             source="google", labels=c("ITS1", "ITS2"),
                             lat.col="Latitude", long.col="Longitude",
                             file=NULL, ext=NULL, height=10, width=12) {
  
  valid.OTU(otu1, otu2)
  .valid.meta(otu1, otu2, meta=meta)
  
  if (is.null(otu2)) {
    num.otus <- 1
  } else {
    num.otus <- 2
  }
  
  .valid.labels(num.otus, labels)
  
  save <- !is.null(file)
  
  if (!source %in% c("google", "osm")) {
    stop("'source' must be one of 'google' or 'osm'.")
  }
  # .valid.factors will check this later, but we do it here to give descriptive
  # error message
  if (!lat.col %in% names(meta)) {
    stop("'lat.col' was not found in meta (see ?sample.locations for help).")
  }
  
  if (!long.col %in% names(meta)) {
    stop("'long.col' was not found in meta (see ?sample.locations for help).")
  }
  
  if (source == "google") {
    if (zoom < 3 || zoom > 21) {
      stop("when source == 'google', zoom must be between 3 and 21 inclusive.")
    }
  } else if (source == "osm") {
    if (zoom < 3 || zoom > 18) {
      stop("when source == 'osm', zoom must be between 3 and 18 inclusive.")
    }
  }
  
  columns <- c(Latitude=lat.col, Longitude=long.col)
  
  # to be used in the ggmap call later
  points.aes <- aes_string(x="Longitude", y="Latitude", size="Counts")
  
  if (!is.null(factor)) {
    
    if (is.null(attr(factor, which="names"))) {
      stop("'factor' must be a named character vector.")
    }
    
    columns <- c(columns, factor)
    
    points.aes <- c(points.aes, aes_string(colour=names(factor)))
    class(points.aes) <- "uneval" # the call to c above strips the class info
  }
  
  if (num.otus == 2) {
    points.aes <- c(points.aes, aes_string(shape="Region"))
    class(points.aes) <- "uneval" # the call to c above strips the class info
  }
  
  meta.data <- .valid.factors(meta, columns, 
                              min.factors=length(columns), max.factors=3)
  
  index <- 1
  for (elem in list(otu1, otu2)) {
    if (is.null(elem)) { break }
    
    meta.data <- cbind(meta.data, colSums(elem[ ,-dim(elem)[2]]))
    # the first two columns are lat/long data, then possibly a metadata factor
    names(meta.data)[index + 2 + !is.null(factor)] <- labels[index]
                       
    index <- index + 1
  }
  
  meta.data <- melt(meta.data, id.vars=c(lat.col, long.col, names(factor)), 
                    variable.name="Region", value.name="Counts")
  
  meta.data.agg <- aggregate(Counts ~ ., data=meta.data, FUN=sum)
  
  buffer <- 0.5
  
  left <- min(meta.data$Longitude) - buffer
  bottom <- min(meta.data$Latitude) - buffer
  right <- max(meta.data$Longitude) + buffer
  top <- max(meta.data$Latitude) + buffer
  
  # create bounding box for map, with a little extra room
  bounding <- c(left=left, bottom=bottom, right=right, top=top)
  
  get_map.args <- list(location=bounding, zoom=zoom, maptype="roadmap", 
                       color="bw", source=source)
  
  if (source == "osm") {
    get_map.args$scale <- OSM_scale_lookup(zoom)
  }
  
  map <- do.call(get_map, get_map.args)
  
  p <- ggmap(map, extent="device") +
       geom_point(points.aes, data=meta.data.agg, alpha=0.7) +
       scale_color_brewer(palette = "Set1") +
       scale_size(range=c(2,10))
  
  if (save) {
    .ggsave.helper(file, ext, width, height, plot=p)
  } else {
    p 
  }
}
