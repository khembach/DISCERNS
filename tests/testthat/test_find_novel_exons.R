context("Predict novel exons")

library(exondiscovery)

gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)

sj <- system.file("extdata", "selected.SJ.out.tab", 
                  package = "exondiscovery", mustWork = TRUE)
bam <- system.file("extdata", "selected.bam", 
                   package = "exondiscovery", mustWork = TRUE)


novel_exon_df <- find_novel_exons(sj_filename = sj, annotation = anno, 
                                  min_unique = 1, bam = bam)

test_that("output from find_novel_exons is correct", {
  expect_s3_class(novel_exon_df, "data.frame")
  expect_equal(ncol(novel_exon_df), 10)
  expect_equal(colnames(novel_exon_df), c("seqnames", "start", "end", "strand",
                                         "lend", "rstart", "unique_left", 
                                         "unique_right", "min_reads", "ID"))
})

test_that("parameters work", {
  ## The parameters should not be <= 0
  expect_error(find_novel_exons(sj_filename = sj, annotation = anno, 
                                min_unique = 1, bam = bam, read_length = -100))
  expect_error(find_novel_exons(sj_filename = sj, annotation = anno, 
                                min_unique = 1, bam = bam, overhang_min = -12))
  expect_error(find_novel_exons(sj_filename = sj, annotation = anno, 
                                min_unique = 1, bam = bam, 
                                min_intron_size = -21))
})