context("Ordination Functions")

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2, meta)

test_that("assist.[cca|rda]", {
  
  assist.functions <- mget(c("assist.cca", "assist.rda"), 
                           envir = as.environment("package:RAM"))
  
  for (func in assist.functions) {
    
    expect_error(func(bad.otu, meta=meta, rank="p"), "valid\\.OTU")
    expect_error(func(ITS1, ITS2, meta=meta[1:10,], rank="p"), "valid\\.meta")
    expect_error(func(ITS1, meta=meta, full=NA, rank="p"))
    expect_error(func(ITS1, meta=meta, exclude=NA, rank="p"))
    expect_error(func(ITS1, meta=meta, rank="x"), "valid.rank")
    
    # commented out as they take too long to run
#     output <- func(ITS1, ITS2, meta, full=TRUE, exclude=2)
#     
#     expect_is(output, "list")
  }
})