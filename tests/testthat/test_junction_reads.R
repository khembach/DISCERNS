context("Predict exons from junction reads and read pairs")

library(exondiscovery)

gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)
bam <- system.file("extdata", "selected.bam", 
                   package = "exondiscovery", mustWork = TRUE)

junc_reads <- filter_junction_reads(bam, lib_type = "PE", stranded = "reverse")
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