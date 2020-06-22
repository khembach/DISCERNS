context("Import novel SJ reads")

library(DISCERNS)

gtf <- system.file("extdata", "selected.gtf", package = "DISCERNS", 
                   mustWork = TRUE)
anno <- prepare_annotation(gtf)
sj <- system.file("extdata", "selected.SJ.out.tab", 
                  package = "DISCERNS", mustWork = TRUE)
bam <- system.file("extdata", "selected.bam", 
                   package = "DISCERNS", mustWork = TRUE)

sj <- fread(sj)
colnames(sj) <- c("seqnames", "start", "end", "strand", "motif", "annotated",
                  "unique", "mutimapping", "maxoverhang")
sj$strand <- c("*", "+", "-")[sj$strand + 1]
sj_gr <- GRanges(sj)
sj_ann <- subsetByOverlaps(sj_gr, anno[["introns"]], type = "equal")
sj_unann <- sj_gr[!(sj_gr %in% sj_ann)]
sj_unann <- sj_unann[(start(sj_unann) - 1) %in% end(anno[["exons"]]) |
                       (end(sj_unann) + 1) %in% start(anno[["exons"]]), ]

reads <- import_novel_sj_reads(bam, sj_unann)

test_that("SJ reads are imported in correct format", {
  expect_s4_class(reads, "GAlignments")
  expect_named(mcols(reads), c("qname", "which_label"))
})