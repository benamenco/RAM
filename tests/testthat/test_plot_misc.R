context("Miscellaneous Plotting Functions")

# set up testing objects
bad.otu <- data.frame(1:10, 1:10, 1:10)
data(ITS1, ITS2, meta)

test_that("group.abundance:input", {
  expect_error(group.abundance(bad.otu, ITS2, rank="p"))
  expect_error(group.abundance(ITS1, ITS2, rank="x"))
  expect_error(group.abundance(ITS1, bad.otu, rank="g"))
  expect_error(group.abundance(ITS1, ITS2, rank="g", top=-1))
  expect_error(group.abundance(ITS1, ITS2, rank="g", top=10,
                               drop.unclassified=NA))
  expect_error(group.abundance(ITS1, ITS2, rank="g", count=NA))
  expect_error(group.abundance(ITS1, ITS2, rank="g", labels="only one label"))
  
  expect_warning(group.abundance(ITS1, ITS2, rank="p", 
                                 labels=c("way", "too", "many")))
})

test_that("group.abundance:output", {
  expect_is(group.abundance(ITS1, ITS2, rank="g", top=10, count=TRUE,
                            drop.unclassified=TRUE, 
                            labels=c("First", "Second"),
                            bw=TRUE, ggplot2=TRUE), 
            "ggplot")
  
  # base graphics
  expect_null(group.abundance(ITS1, ITS2, rank="g", top=10, count=FALSE,
                            drop.unclassified=TRUE, labels=c("OTU1", "OTU2"),
                            bw=FALSE, ggplot2=FALSE))
  dev.off()
})

test_that("group.indicators:input", {
  expect_error(group.indicators(bad.otu, ITS2, meta=meta, factor=c(Crop="Crop"), 
                                rank="p"))
  expect_error(group.indicators(ITS1, ITS2, meta=meta, factor=c(Crop="Crop"), rank="x"))
  expect_error(group.indicators(ITS1, ITS2, bad.otu, meta=meta, factor=c(Crop="Crop"),
                                rank="g"))
  
  expect_error(group.indicators(ITS1, ITS2, meta=meta, factor=c(Crop="Crop"), 
                                thresholds=c(0.1, 0.1, 0.1, -0.1), rank="g"))
  
  expect_error(group.indicators(ITS1, ITS2, meta=meta, factor=c(Crop="Crop"),
                                labels=c("not enough"), rank="g"))
  
  expect_error(group.indicators(ITS1, ITS2, meta=meta, factor=c(Crop="Croz"), rank="g"))

  expect_warning(group.indicators(ITS1, meta=meta, factor=c(Crop="Crop"),
                                  labels=c("TOO", "MANY"), rank="g"))
})

test_that("group.indicators:output", {
  
  expect_is(group.indicators(ITS1, ITS2, meta=meta, factor=c(Prov.="Province"),
                             rank="g", thresholds=c(0.9, 0.7, 0.5, 0.05),
                             labels=c("First OTU", "Second OTU")),
            "ggplot")
})

test_that("group.spatial:input", {
  
  expect_error(group.spatial(ITS1))
  expect_error(group.spatial(ITS1, meta=meta, date.col="BAD_DATE", 
                             province.col="Province",
                             rank="p", group="Ascomycota"),
               "date.col")
  
  expect_error(group.spatial(bad.otu, meta=meta, date.col="Harvestdate", 
                             province.col = "Province",
                             rank="p", group="Ascomycota"),
               "valid\\.OTU")
  
  expect_error(group.spatial(ITS1, meta=meta, date.col="Harvestdate", 
                             province.col="Province", rank="p",
                             group="Ascomycota",
                             breaks="BAD_BREAKS"),
               "breaks")
  
  expect_error(group.spatial(ITS1, meta=meta, date.col="Harvestdate",
                             rank="p", group="Ascomycota",
                             province.col="BAD_PROVINCE"),
               "province")
  
  expect_error(group.spatial(ITS1, meta=meta, date.col="Harvestdate", 
                             rank="x", group="Ascomycota",
                             province.col="Province"),
               "\\.valid\\.rank")
  
  expect_error(group.spatial(ITS1, meta=meta, date.col="Harvestdate",
                             rank="p", group="BAD_GROUP", 
                             province.col = "Province"),
               "\\.get\\.tax\\.group")
})

# this case is commented out due to time restrictions
# test_that("group.spatial:output", {
#   
#   expect_is(group.spatial(ITS1, meta=meta, date.col="Harvestdate", 
#                           province.col="Province", rank="p", 
#                           group=c("Ascomycota", "Basidiomycota"),
#                           breaks=2),
#             "ggplot")
# })

test_that("group.temporal:input", {
  
  expect_error(group.temporal(ITS1))
  expect_error(group.temporal(ITS1, meta=meta, date.col="BAD_DATE", 
                             rank="p", group="Ascomycota"),
               "date.col")
  
  expect_error(group.temporal(bad.otu, meta=meta, date.col="Harvestdate", 
                             rank="p", group="Ascomycota"),
               "valid\\.OTU")
  
  expect_error(group.temporal(ITS1, meta=meta, date.col="Harvestdate",
                             rank="p", group="Ascomycota",
                             factors=c(Ergo="BAD_NAME")),
               "\\.valid\\.factors")
  
  expect_error(group.temporal(ITS1, meta=meta, date.col="Harvestdate",
                             rank="p", group="BAD_GROUP",
                             factors=c(Ergosterol="Ergosterol_ppm")),
               "\\.get\\.tax\\.group")
  
  expect_error(group.temporal(ITS1, meta=meta, date.col="Harvestdate",
                              rank="x", group="Ascomycota",
                              factors=c(Ergosterol="Ergosterol_ppm")),
               "\\.valid\\.rank")
})

test_that("group.temporal:output", {
  
  # this isn't a ggplot object because we're manually arranging 2 plots on 1 grid
  expect_null(group.temporal(ITS1, meta=meta, date.col="Harvestdate",
                           rank="p", group="Ascomycota",
                           factors=c(Ergosterol="Ergosterol_ppm")))
  dev.off()
})

test_that("group.top.samples:input", {
  for (func in mget(c("group.top.number", "group.top.percent"), 
                    envir=as.environment("package:RAM"))) {
    
    expect_error(func(bad.otu, ITS2), "valid\\.OTU")
    expect_error(func(ITS1, ITS2, top=-1), "top")
    expect_error(func(ITS1, ITS2, top="baaaaad"), "top")
    expect_error(func(ITS1, ITS2, drop.unclassified = NA))
    expect_error(func(ITS1, ITS2, labels="TOO_FEW"), "labels")
    expect_warning(func(ITS1, top=3, labels=c("TOO", "MANY")), "labels")
    expect_error(func(ITS1, ITS2, ggplot2=NA))
    expect_error(func(ITS1, ITS2, bw=NA))
  }
})

test_that("group.top.number:input", {
  
  expect_warning(group.top.number(ITS1, top=100000))
})

test_that("group.top.number:output", {
  
  expect_is(group.top.number(ITS1, ITS2, top=5, drop.unclassified=TRUE,
                             labels=c("test1", "test2"), bw=FALSE, ggplot2=TRUE),
            "ggplot")
  
  expect_null(group.top.number(ITS1, top=3, drop.unclassified=FALSE,
                               labels="Hello", bw=TRUE, ggplot2=FALSE))
  dev.off()
})

test_that("group.top.percent:input", {
  
  expect_error(group.top.percent(ITS1, top=101))
})

test_that("group.top.percent:output", {
  
  expect_is(group.top.number(ITS1, top=2, drop.unclassified=FALSE,
                             labels="the only one", bw=TRUE, ggplot2=TRUE),
            "ggplot")
  
  expect_null(group.top.number(ITS1, ITS2, top=4, drop.unclassified=TRUE,
                               labels=c("gg", "nore"), bw=FALSE, ggplot2=FALSE))
  dev.off()
})

test_that("pcoa.plot:input", {
  
  expect_error(pcoa.plot(bad.otu, meta, c(Crop="Crop"), rank="p"),
               "valid\\.OTU")
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="BAD_RANK"),
               "valid\\.rank")
  expect_error(pcoa.plot(ITS1, meta[1:10, ], c(Crop="Crop"), rank="p"),
               "valid\\.meta")
  expect_error(pcoa.plot(ITS1, meta, c(bad="BAD_FACTOR"), rank="p"),
               "valid\\.factors")
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", ggplot2=NA))
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", 
                         sample.labels=NA))
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", ggplot2=FALSE,
                         ellipse=2),
               "ellipse")
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p",
                         ggplot2=FALSE, ellipse=NA))
  
  expect_warning(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", ggplot2=TRUE,
                         ellipse=1))
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", 
                         dist.method="BAD_DIST_METHOD"))
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p",
                         stand.method="BAD_STAND_METHOD"))
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p",
                         stand.method=NA))
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", top=-1),
               "top")
  
  expect_error(pcoa.plot(ITS1, meta, c(Crop="Crop"), rank="p", bw=NA))
})

test_that("pcoa.plot:output", {
  
  expect_is(pcoa.plot(ITS1, meta, rank="p", factors=c(Crop="Crop")),
            "ggplot")
  
  expect_null(pcoa.plot(ITS2, meta, rank="p", factors=c(City="City"), top=0,
                        ggplot2=FALSE, bw=TRUE, ellipse=1, sample.labels=FALSE))
})