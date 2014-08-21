context("Taxonomy Functions")

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2)

test_that("get.rank:input", {
     
 expect_error(get.rank(bad.otu, rank="c"), "valid\\.OTU")
 expect_error(get.rank(ITS1, bad.otu, rank="c"), "valid\\.OTU")
 expect_error(get.rank(ITS1, rank=NA), "valid\\.rank")
 expect_error(get.rank(ITS1, rank=c("k", "p", "c")), "valid\\.rank")
 
 expect_error(get.rank(ITS1, rank="x"), "valid\\.rank")
 expect_error(get.rank(bad.otu, bad.otu, rank="x"), "valid\\.OTU|valid\\.rank")
})

test_that("get.rank:output", {
  
  single.output <- get.rank(ITS1, rank="c")
  full.output <- get.rank(ITS1, ITS2)
  
  # check output data type
  expect_is(single.output, "data.frame")
  expect_is(full.output, "list")
  
  # check output data is valid
  expect_is(valid.OTU(single.output), "NULL")
  expect_is(valid.OTU(full.output$otu2$genus), "NULL")
  # check output data dimensions
  expect_equal(dim(single.output)[2], dim(ITS1)[2])
  expect_equal(length(full.output), 2)
  expect_equal(length(full.output$otu2), 7)
  # check all columns (other than last) are numeric
  tax <- dim(single.output)[2]
  expect_true(all(sapply(single.output[ ,-tax], FUN=is.numeric)))
  
  expect_identical(single.output, full.output$otu1$class)
  })

test_that("tax.split:input", {
  
  expect_error(tax.split(bad.otu, rank="c"), "valid\\.OTU")
  expect_error(tax.split(ITS1, bad.otu, rank="c"), "valid\\.OTU")
  expect_error(tax.split(ITS1, rank=NA), "valid\\.rank")
  expect_error(tax.split(ITS1, rank=c("k", "p", "c")), "valid\\.rank")
  
  expect_error(tax.split(ITS1, rank="x"), "valid\\.rank")
  expect_error(tax.split(bad.otu, bad.otu, rank="x"), "valid\\.OTU|valid\\.rank")
})

test_that("tax.split:output", {
  
  single.output <- tax.split(ITS1, rank="c")
  full.output <- tax.split(ITS1, ITS2)
  
  # check output data type
  expect_is(single.output, "data.frame")
  expect_is(full.output, "list")
  
  # check output data is valid
  expect_equal(dim(single.output), dim(ITS1))
  # check output data dimensions
  expect_equal(length(full.output), 2)
  expect_equal(length(full.output$otu2), 7)
  # check all columns (other than last) are numeric
  tax <- dim(single.output)[2]
  expect_true(all(sapply(single.output[ ,-tax], FUN=is.numeric)))
  
  expect_identical(single.output, full.output$otu1$class)
})

test_that("tax.abund:input", {
  
  expect_error(tax.abund(bad.otu, rank="c"), "valid\\.OTU")
  expect_error(tax.abund(ITS1, bad.otu, rank="c"), "valid\\.OTU")
  expect_error(tax.abund(ITS1, rank=NA), "valid\\.rank")
  expect_error(tax.abund(ITS1, rank=c("k", "p", "c")), "valid\\.rank")
  
  expect_error(tax.abund(ITS1, rank="x"), "valid\\.rank")
  expect_error(tax.abund(bad.otu, bad.otu, rank="x"), "valid\\.OTU|valid\\.rank")
})

test_that("tax.abund:output", {
  
  single.output <- tax.abund(ITS1, rank="c")
  full.output <- tax.abund(ITS1, ITS2)
  
  # check output data type
  expect_is(single.output, "data.frame")
  expect_is(full.output, "list")
  # check output data dimensions
  expect_equal(length(full.output), 2)
  expect_equal(length(full.output$otu2), 7)
  # check all columns are numeric
  expect_true(all(sapply(single.output, FUN=is.numeric)))
  
  expect_identical(single.output, full.output$otu1$class)
})

test_that("tax.fill:input", {
  
  expect_error(tax.fill(ITS1, NA), "downstream")
  expect_error(tax.fill(bad.otu, TRUE), "valid\\.OTU")
  expect_error(tax.fill(ITS1, c(TRUE, FALSE, TRUE)), "downstream")
  expect_error(tax.fill(ITS2, "NA"), "downstream")
})

test_that("tax.fill:output", {
  
  output <- tax.fill(ITS1, TRUE)
  
  expect_is(output, "data.frame")
  expect_is(output$taxonomy, "character")
  
  expect_equal(dim(output), dim(ITS1))
  # ensure that output passes valid.OTU check, which returns NULL if OTU is good
  expect_null(valid.OTU(output))
})