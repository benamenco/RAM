.true.calc <- function(otu1, otu2=NULL, index, mode) {
  given.both <- !is.null(otu2) 
  
  valid.OTU(otu1, otu2)
  otu1.t <- transpose.OTU(otu1)
  
  if (!is.character(index)) {
    stop("argument 'index' must be either 'simpson' or 'shannon'.")
  }
  if (length(index) > 1) {
    warning("only the first element of argument 'index' is being used.")
    index <- index[1]
  }
  
  if (given.both) { 
    otu2.t <- transpose.OTU(otu2)
  }
  
  if (mode == "diversity") {
    func <- .true.div
  } else if (mode == "evenness") {
    func <- .true.even
  }
  
  methods <- c("simpson", "shannon")
  ### is this error message descriptive enough?
  index <- match.arg(index, methods)
  
  # get the sample names; remove the last name (it will be taxonomy)
  cnames <- names(otu1)[-length(names(otu1))]
  
  if (given.both) {
    res1 <- func(otu1.t, index)
    res2 <- func(otu2.t, index)
    
    rnames <- c("otu1", "otu2")
    
    output <- matrix(rbind(res1, res2), nrow=2, dimnames=list(rnames, cnames))
    return(output)
    
  } else {
    res <- func(otu1.t, index)
    
    rnames <- c("otu1")
    matrix(res, nrow=1, dimnames=list(rnames, cnames))
  }
}

true.diversity <- function(otu1, otu2=NULL, index="simpson") {
  .true.calc(otu1, otu2, index=index, mode="diversity")
}

evenness <- function(otu1, otu2=NULL, index="simpson") {
  .true.calc(otu1, otu2, index=index, mode="evenness")
}

.true.div <- function(OTU, index) {
  
  if (index == "simpson") {
    return(diversity(OTU, index="invsimpson"))
    
  } else if (index == "shannon") {
    return(exp(diversity(OTU, index="shannon", base = exp(1))))
  }
}

.true.even <- function(OTU, index) {
  
  if (index == "simpson") {
    return(diversity(OTU, index="invsimpson") / specnumber(OTU))
    
  } else if (index == "shannon") {
    return(diversity(OTU, index="shannon", base = exp(1)) / 
              log(specnumber(OTU), base=exp(1)))
  }
}