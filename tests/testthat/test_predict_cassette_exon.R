context("Predict cassette exons")

library(exondiscovery)

gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)

## Two SJs that define a cassette exon
sj_unann <- GRanges("19", IRanges(c(10, 17119370, 17120131), c(20, 17120112, 17122273)), 
                    strand = "-")
ce <- predict_cassette_exon(sj_unann, anno[["introns"]])

test_that("output from predict_cassette_exon is correct", {
  expect_type(ce, "list")
  expect_named(ce, c("ne", "sj"))
  expect_s3_class(ce[["ne"]], "data.frame")
  expect_s4_class(ce[["sj"]], "GRanges")

  ## coordinates of the novel cassette exon
  expect_equal(ce[["ne"]]$lend, 17119369)
  expect_equal(ce[["ne"]]$start, 17120113)
  expect_equal(ce[["ne"]]$end, 17120130)
  expect_equal(ce[["ne"]]$rstart, 17122274)
  
  ## single SJ as input --> no cassette exon
  ce <- predict_cassette_exon(sj_unann[1], anno[["introns"]])
  expect_equal(nrow(ce[["ne"]]), 0)
  expect_length(ce[["sj"]], 1)
})

