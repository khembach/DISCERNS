context("Extend GTF annotation")

library(exondiscovery)

gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery", 
                   mustWork = TRUE)

novel_exons <- data.frame(seqnames = c("19", "22"),
                          start = c(47228064,41737092),
                          end = c(47228185, 41737150), strand = c("-", "+"),
                          lend = c(47226541, 41736141),
                          rstart = c(47228589, 41738533), ID = c(1, 2))
gtf <- import(gtf)
new_gtf <- extend_gtf(gtf, novel_exons)

test_that("extend_gtf() creates correct output",{
  expect_s4_class(new_gtf, "GRanges")
  expect_equal(names(mcols(new_gtf)), names(mcols(gtf)))
  ## Check that we have at least one novel exon entry
  expect_true(sum(grepl("_2", new_gtf$exon_id)) > 0)
  ## Check that the two novel exons are not in gtf
  expect_length(subsetByOverlaps(unique(gtf), 
                                 GRanges(novel_exons), type = "equal"), 0)
  ## Check that the two novel exons have been added to new_gtf
  expect_length(subsetByOverlaps(unique(new_gtf), 
                                 GRanges(novel_exons), type = "equal"), 2)
  ## Check that the novel exon entries have the novel ID appended to all
  ## transcript IDa and names
  novel_entries <- subsetByOverlaps(new_gtf, 
                                    GRanges(novel_exons)[1], type = "equal")
  expect_equal(sum(grepl("_1", novel_entries$transcript_id)), 
               length(novel_entries))
  expect_equal(sum(grepl("_1", novel_entries$transcript_name)), 
               length(novel_entries))
  expect_equal(sum(grepl("_1", novel_entries$exon_id)), 
               length(novel_entries))
  ## Check that exon_number and exon_version are set to NA
  expect_true(all(is.na(novel_entries$exon_number)))
  expect_true(all(is.na(novel_entries$exon_version)))
})
  