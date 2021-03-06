context("Prepare GTF annotation")

library(DISCERNS)

gtf <- system.file("extdata", "selected.gtf", package = "DISCERNS", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)

test_that("annotation is prepared correctly", {
  expect_type(anno, "list")
  expect_named(anno, c("exons", "introns", "txdb"))
  expect_is(anno[["exons"]], "GRanges")
  expect_is(anno[["introns"]], "GRanges")
  expect_is(anno[["txdb"]], "TxDb")
  
  gtf_file <- ""
  expect_error(prepare_annotation(""), paste0("File ", gtf_file, 
                                              " does not exist."))
})
