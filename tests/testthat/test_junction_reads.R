context("Predict exons from junction reads and read pairs")

library(DISCERNS)

gtf <- system.file("extdata", "selected.gtf", package = "DISCERNS", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)
bam <- system.file("extdata", "selected.bam", 
                   package = "DISCERNS", mustWork = TRUE)

junc_reads <- filter_junction_reads(bam, lib_type = "PE", stranded = "reverse")

test_that("filter_junction_reads works correctly", {
  expect_error(filter_junction_reads(bam, lib_type = "pe"), 
               'Parameter lib_type has to be either "SE" or "PE".')
  expect_error(filter_junction_reads(bam, stranded = "not stranded"), 
               'Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  
  expect_type(junc_reads, "list")
  expect_true(length(junc_reads) > 0)
})

junc_reads <- do.call(c, junc_reads) 
read_pred <- predict_jr_exon(junc_reads, anno)

test_that("predict_jr_exon() input and output are correct", {
  expect_error(predict_jr_exon(junc_reads = GAlignments(), annotation = anno))
  
  expect_s3_class(read_pred, "data.frame")
})


read_pred <- predict_jrp_exon(junc_reads, anno)

test_that("predict_jrp_exon() input, output and parameters work", {
  expect_error(predict_jrp_exon(junc_reads = GAlignments(), annotation = anno))

  expect_s3_class(read_pred, "data.frame")
  
  ## Some parameters must not be <= 0
  expect_error(predict_jrp_exon(junc_reads = junc_reads, annotation = anno, 
                                min_unique = 1, bam = bam, read_length = -100))
  expect_error(predict_jrp_exon(junc_reads = junc_reads, annotation = anno, 
                                min_unique = 1, bam = bam, overhang_min = -12))
  expect_error(predict_jrp_exon(junc_reads = junc_reads, annotation = anno, 
                                min_unique = 1, bam = bam, 
                                min_intron_size = -21))
})