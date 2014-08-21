.analysis.helper <- function(otu1, otu2=NULL, meta, full, exclude, mode, rank) {
  valid.OTU(otu1, otu2)
  .valid.meta(otu1, otu2, meta)
  .valid.rank(rank)
  
  single.otu <- is.null(otu2)
  
  if (!single.otu) {
    return(list(otu1=.analysis.helper(otu1, meta=meta, full=full,
                                      exclude=exclude, mode=mode, rank=rank),
                otu2=.analysis.helper(otu2, meta=meta, full=full,
                                      exclude=exclude, mode=mode, rank=rank)))
  }
  
  tax<-tax.abund(otu1, otu2=NULL, rank=rank, drop.unclassified=TRUE,
                 top=NULL, count=TRUE, mode="number")
  
  
  if (is.null(exclude)) {
    meta.remain <- meta
  } else {
    meta.remain <- meta[ , -exclude]
  }
  
  if (full) {
    form <- formula(tax ~ .)
  } else {
    form <- formula(tax ~ 1)
  }
  
  if (mode == "cca") {
    mod <- vegan::cca(form, data=meta.remain)
  } else if (mode == "rda") {
    mod <- vegan::rda(form, data=meta.remain)
  }
  
  good <- goodness(mod,summ = TRUE)
  vif.scores <- vif.cca(mod)
  pct.var <- summary(mod)$concont$importance
  CCA.eig <- mod$CCA$eig
  CA.eig <- mod$CA$eig
  anv <- vegan::anova.cca(mod)

  return(list(GOF=good, VIF=vif.scores, percent_variation=pct.var, CCA_eig=CCA.eig, 
              CA_eig=CA.eig, anova=anv))
}

assist.cca <- function(otu1, otu2=NULL, meta, full=TRUE, exclude=NULL, rank) {
  .analysis.helper(otu1, otu2, meta, full, exclude, mode="cca", rank=rank)
}


assist.rda <- function(otu1, otu2=NULL, meta, full=TRUE, exclude=NULL, rank) {
  .analysis.helper(otu1, otu2, meta, full, exclude, mode="rda", rank=rank)
}
