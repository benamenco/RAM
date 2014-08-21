context("Miscellaneous Functions")

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2, meta)

test_that("percent.classified", {
  
  expect_error(percent.classified(bad.otu, "p"), "valid\\.OTU")
  expect_error(percent.classified(ITS1, "x"), "valid\\.rank")
  expect_is(percent.classified(ITS2, "c"), "numeric")
  expect_warning(percent.classified(ITS1, "s"), "no OTUs")
})

test_that("diversity.indices", {
  
  diversity.functions <- mget(c("true.diversity", "evenness"),
                              envir = as.environment("package:RAM"))
  
  for (func in diversity.functions) {
    
    expect_error(func(bad.otu), "valid\\.OTU")
    expect_error(func(ITS1, bad.otu), "valid\\.OTU")
    expect_error(func(ITS1, index="BAD_INDEX"))
    
    output <- func(ITS1, ITS2, index="simpson")
    
    expect_is(output, "matrix")
    expect_true(is.numeric(output))
    expect_equal(dim(output), c(2, dim(ITS1)[2] - 1))
  }
})

test_that("shared.OTU", {
  
  expect_error(shared.OTU(bad.otu), "valid\\.OTU")
  
  output <- shared.OTU(ITS1)
  
  expect_is(output, "list")
  expect_equal(length(output), 8L)
})