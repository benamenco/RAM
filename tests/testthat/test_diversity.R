context("Diversity Functions")

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2)

test_that("true.diversity:input", {
  
  expect_error(true.diversity(bad.otu, index="simpson"), "valid\\.OTU")
  expect_error(true.diversity(ITS1, bad.otu, index="simpson"), "valid\\.OTU")
  expect_error(true.diversity(ITS1, index=NA), "\\.true\\.calc")
  expect_error(true.diversity(ITS1, index=bad.otu), "\\.true\\.calc")
  expect_error(true.diversity(ITS1, index="s"), "match\\.arg")
  expect_warning(true.diversity(ITS1, index=c("simp", "shan")), "only the first element")
})

test_that("true.diversity:output", {
  
  single.output <- true.diversity(ITS1, index="simpson")
  double.output <- true.diversity(ITS1, ITS2, index="shannon")
  num.samples <- dim(ITS1)[2] - 1
  
  expect_is(single.output, "matrix")
  expect_is(double.output, "matrix")
  
  expect_equal(dim(single.output), c(1, num.samples))
  expect_equal(dim(double.output), c(2, num.samples))
  
  expect_equal(colnames(single.output), names(ITS1)[1:num.samples])
  expect_equal(colnames(double.output), names(ITS1)[1:num.samples])
  
  expect_equal(rownames(single.output), c("otu1"))
  expect_equal(rownames(double.output), c("otu1", "otu2"))
})

test_that("evenness:input", {
  
  expect_error(evenness(bad.otu, index="simpson"), "valid\\.OTU")
  expect_error(evenness(ITS1, bad.otu, index="simpson"), "valid\\.OTU")
  expect_error(evenness(ITS1, index=NA), "\\.true\\.calc")
  expect_error(evenness(ITS1, index=bad.otu), "\\.true\\.calc")
  expect_error(evenness(ITS1, index="s"), "match\\.arg")
  expect_warning(evenness(ITS1, index=c("simp", "shan")), "only the first element")
})

test_that("evenness:output", {
  
  single.output <- evenness(ITS1, index="simpson")
  double.output <- evenness(ITS1, ITS2, index="shannon")
  num.samples <- dim(ITS1)[2] - 1
  
  expect_is(single.output, "matrix")
  expect_is(double.output, "matrix")
  
  expect_equal(dim(single.output), c(1, num.samples))
  expect_equal(dim(double.output), c(2, num.samples))
  
  expect_equal(colnames(single.output), names(ITS1)[1:num.samples])
  expect_equal(colnames(double.output), names(ITS1)[1:num.samples])
  
  expect_equal(rownames(single.output), c("otu1"))
  expect_equal(rownames(double.output), c("otu1", "otu2"))
})