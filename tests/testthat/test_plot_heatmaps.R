context("Heatmap Plotting Functions")

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2, meta)

test_that("group.heatmap:input", {
  
  expect_error(group.heatmap(bad.otu, meta=meta, rank="g"), "valid\\.OTU")
  expect_error(group.heatmap(ITS1, meta=meta, factors=c(Crop="Crop", Ergo=1),
                             rank="g"), 
               "\\.valid\\.factors")
  expect_error(group.heatmap(ITS1, meta[1:10, ], factors=c(Crop="Crop"),
                             rank="g"),
               "\\.valid\\.meta")
  
  expect_error(group.heatmap(ITS1, meta=meta, factors=c(Crop="Crop"),
                             rank="g", drop.unclassified=NA))
  
  expect_error(group.heatmap(ITS1, meta=meta, factors=c(Crop="Crop"),
                             rank="g", count=NA))
  
  expect_error(group.heatmap(ITS1, meta=meta, factors=c(Crop="Crop"),
                             rank="g", top=-1))
  
  expect_warning(group.heatmap(ITS1, meta=meta, factors=c(Crop="Crop"),
                               rank="g", cut=100000000))
  dev.off()
  
  expect_warning(group.heatmap(ITS1, meta=meta, factors="Crop", rank="g"))
  dev.off()
  
})

test_that("group.heatmap:output", {
  
  # make sure a valid call works
  expect_is(group.heatmap(ITS1, meta=meta, rank="g", top=10, 
                          factors=c(Crop="Crop", Ergo="Ergosterol_ppm"),
                          drop.unclassified=TRUE, cut=0.5, count=FALSE), 
            "annHeatmap")
  dev.off()
})

test_that("group.heatmap.simple:input", {
  
  expect_error(group.heatmap.simple(bad.otu, meta=meta, rank="g"), "valid\\.OTU")
  expect_error(group.heatmap.simple(ITS1, meta=meta, row.factor=c(Crop="Crozz"),
                             rank="g"), 
               "\\.valid\\.factors")
  
  expect_error(group.heatmap.simple(ITS1, meta[1:10, ], rank="g"),
               "\\.valid\\.meta")
  
  expect_error(group.heatmap.simple(ITS1, meta=meta, rank="g", top=-1))
  
  expect_error(group.heatmap.simple(ITS1, meta=meta, rank="g", top=10,
                                    count=NA))
  
  expect_error(group.heatmap.simple(ITS1, meta=meta, rank="g", top=10,
                                    drop.unclassified=NA))
  
  expect_error(group.heatmap.simple(ITS1, meta=meta, rank="g", top=10,
                                    dendro=TRUE))
  
  expect_warning(group.heatmap.simple(ITS1, meta=meta, rank="g", top=10,
                                      row.factor=c(Prov="Province"), dendro="both"),
                 "clustering will not occur")
  dev.off()
})

test_that("group.heatmap.simple:output", {
  
  expect_null(group.heatmap.simple(ITS1, meta=meta, rank="o", top=7,
                                   row.factor=c(Crop="Crop"),
                                   drop.unclassified = FALSE,
                                   dendro="none"))
  dev.off()
})

test_that("dissim.heatmap:input", {
  
  expect_error(dissim.heatmap(bad.otu, meta=meta), "valid\\.OTU")
  expect_error(dissim.heatmap(ITS1, meta=meta, row.factor=c(Crop="Crozz")), 
               "\\.valid\\.factors")
  
  expect_error(dissim.heatmap(ITS1, meta=meta[1:10, ], row.factor=c(Crop="Crop")),
               "\\.valid\\.meta")
  
  expect_error(dissim.heatmap(ITS1, meta=meta, stand.method=NA))
  
  expect_error(dissim.heatmap(ITS1, meta=meta, dissim.method=NA))
  
  expect_error(dissim.heatmap(ITS1, meta=meta, col.factor=c(1:5)))
  
  expect_warning(dissim.heatmap(ITS1, meta=meta, 
                               row.factor=c(Ergo="Ergosterol_ppm")))
  dev.off()
  
})

test_that("dissim.heatmap:output", {
  
  expect_null(dissim.heatmap(ITS1, meta=meta, row.factor=c(Crop="Crop"),
                             col.factor=c(City="City"), 
                             stand.method="chi.square",
                             dissim.method="euclidean"))
  dev.off()
})