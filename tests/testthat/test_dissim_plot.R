context("Dissimilarity Plots")

# many of the functions in dissim.X.plot take the same parameters, and so we
# apply the same tests to those functions

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2)

# Type 1 functions have the form 
# <funcion>(otu1, otu2=NULL, file=NULL, dist.methods=NULL)

type1 <- list("dissim.clust.plot", "dissim.eig.plot", "dissim.alleig.plot",
              "dissim.GOF.plot", "dissim.tree.plot", "dissim.pvar.plot")

test_that("dissim.type1:input", {
  
  for (name in type1) {
    
    func <- match.fun(name)
    
    expect_error(func(bad.otu), "valid\\.OTU")
    expect_error(func(ITS1, bad.otu), "valid\\.OTU")
    expect_error(func(ITS1, dist.methods=NA), "vegdist")
    expect_error(func(ITS1, dist.methods=TRUE), "vegdist")
    expect_error(func(ITS1, dist.methods=list("spiderman", "morisita")), "vegdist")
    expect_error(func(ITS1, dist.methods=list("morisita", "spiderman")), "vegdist")
    expect_error(func(ITS1, dist.methods=c("morisita", "spiderman")), "vegdist")
    expect_error(func(ITS1, file=""), "\\.ensure\\.filepath")
    expect_error(func(ITS1, file=NA), "\\.ensure\\.filepath")
    expect_error(func(ITS1, file=TRUE), "\\.ensure\\.filepath")
    expect_error(func(ITS1, file=c("hi", "mom")), "\\.ensure\\.filepath")
    expect_error(func(ITS1, file=list("hello", 2)), "\\.ensure\\.filepath")
    expect_error(func(ITS1, file=c("songbird", "morisita")), "\\.ensure\\.filepath")
  }
})

test_that("dissim.others:input", {
  
  clust.funcs <- list("dissim.clust.plot", "dissim.tree.plot")
  
  for (name in clust.funcs) {
    
    func <- match.fun(name)
    
    # this test fails due to hclust not safely handling NA for argument 'method'
    # expect_error(func(ITS1, clust.methods=NA), "hclust")
    expect_error(func(ITS1, clust.methods=TRUE), "hclust")
    expect_error(func(ITS1, clust.methods=list(2, "ward")), "hclust")
  }
  
  ord.funcs <- list("dissim.ord.plot")
  
  for (name in ord.funcs) {
    
    func <- match.fun(name)
    
    expect_error(func(ITS1, dist.methods=NA), "vegdist")
    expect_error(func(ITS1, dist.methods=TRUE), "vegdist")
    expect_error(func(ITS1, dist.methods=list(2, "ward")), "vegdist")
    # this test depends on the dimensions of the test data
    expect_error(func(ITS1, k=17), "specified value for k was")
    expect_error(func(ITS1, k=c("hello", "johnny")), "k must be a numeric")
  }
})

dev.off()