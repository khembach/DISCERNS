library(exondiscovery)

gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)

test_that("annotation is prepared correctly", {
  expect_type(anno, "list")
  expect_named(anno, c("exons", "introns", "txdb"))
  expect_is(anno[["exons"]], "GRanges")
  expect_is(anno[["introns"]], "GRanges")
  expect_is(anno[["txdb"]], "TxDb")
  })
