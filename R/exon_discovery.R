#' Predict novel cassette exons from SJ.out.tab file from STAR
#'
#' @param sj_filename Path to SJ.out.tab file.
#' @param annotation List with exon and intron annotation as GRanges. Created
#'   with prepare_annotation().
#' @param min_unique Minimal number of reads required to map over a splice
#'   junction.
#' @param verbose Print messages about the progress.
#' @param gzipped A logical scalar. Is the input file gzipped? Default FALSE.
#' @param bam Path to BAM file.
#'
#' @return Data frame with novel exons.
#'
#' @importFrom data.table fread
#' @import GenomicRanges
#' @import IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicFeatures genes exonsBy
#' @importFrom dplyr full_join mutate
#'
#' @export
#'
find_novel_exons <- function(sj_filename, annotation, min_unique = 1,
                                      gzipped = FALSE, verbose = TRUE, bam) {
  exons <- annotation[["exons"]]
  introns <- annotation[["introns"]]

  if (verbose) message("Step 1: Reading splice junctions from file")

  if (gzipped) {
    sj <- fread(paste0("zcat ", sj_filename))
  } else {
    sj <- fread(sj_filename)
  }
  colnames(sj) <- c("seqnames", "start", "end", "strand", "motif", "annotated",
                    "unique", "mutimapping", "maxoverhang")

  #strand: (0: undefined, 1: +, 2: -)
  sj$strand <- c("*", "+", "-")[sj$strand + 1]

  ### TODO: create functions for the different prediction "types"
  ## 1) casette exons 2) predictions from single novel SJ
  ## 3) predictions from reads or read pairs with 2 SJ
  ## let the user decide which ones to run

  ## ============= junction prefiltering ============= #

  if (verbose) message("Step 2: Prefiltering splice junctions")

  sj_gr <- GRanges(sj)

  ## filter out all annotated junctions
  sj_ann <- subsetByOverlaps(sj_gr, introns, type = "equal")
  sj_unann <- sj_gr[!(sj_gr %in% sj_ann)]

  ## filter out all junctions with less than min_unique reads
  sj_unann <- sj_unann[mcols(sj_unann)$unique >= min_unique]

  ## filter out all novel junctions that do not touch an annotated exon on the
  ## same strand
  sj_unann <- sj_unann[(start(sj_unann) - 1) %in% end(exons) |
                         (end(sj_unann) + 1) %in% start(exons), ]

  ## ======== Cassette exon prediction ======== #

  if (verbose) message("Step 3: Predicting cassette exons")

  ce <- predict_cassette_exon(sj_unann, introns)
  novel_exons <- ce[["ne"]]
  sj_unann <- ce[["sj"]]

  ## ==== Find novel exon coordinates from single novel splice junctions ======
  if (verbose)
    message("Step 4: Find exon coordinates of exons adjacent to novel splice
            junctions")

  reads <- import_novel_sj_reads(bam, sj_unann)

  ebyTr <- exonsBy(annotation[["txdb"]], by = "tx", use.names = TRUE)
  gtxdb <- genes(annotation[["txdb"]])

  ## convert GRanges to data.frame otherwise we cannot use apply functions
  sj_unann <- data.frame(sj_unann, stringsAsFactors = TRUE)

  ## annotate which end of a splice junction is touching an annotated exon
  touching <- mapply(sj_touching_exon, sj_unann$seqnames, sj_unann$start,
                     sj_unann$end, sj_unann$strand,
                     MoreArgs = list(exons = exons))
  sj_unann$touching <- touching
  ## remove the ambiguous sj that touch different genes on both ends
  sj_unann <- sj_unann[!is.na(touching), ]

  res <- identify_exon_from_sj(sj_unann, reads = reads,
                                txdb = annotation[["txdb"]],
                                gtxdb = gtxdb, ebyTr = ebyTr)

  if (nrow(res) > 0) {
    novel_exons <- rbind(novel_exons, res)
  }

  ## ============ Predict novel exons from reads with 2 junctions ==============
  if (verbose) message("Step 5: Predict novel exons from reads with 2 junctions")

  junc_reads <- filter_junction_reads(bam)
  read_pred <- predict_jr_exon(junc_reads, annotation)

  ## make sure the factors have the same levels
  combined <- union(levels(read_pred$strand), levels(novel_exons$strand))
  ## predictions from both the reads and the SJ.out.tab
  novel_exons <- full_join(mutate(read_pred,
                                  strand = factor(strand, levels = combined)),
                           mutate(novel_exons,
                                  strand = factor(strand, levels = combined)))

  ## ======= Predict novel exons from reads pairs with each one junction =======
  if (verbose) message("Step 6: Predict novel exons from reads pairs with each 1
                      junction")

  read_pair_pred <- predict_jrp_exon(junc_reads, annotation)
  combined <- union(levels(read_pair_pred$strand), levels(novel_exons$strand))
  novel_exons <- full_join(mutate(read_pair_pred,
                                  strand = factor(strand, levels = combined)),
                           mutate(novel_exons,
                                  strand = factor(strand, levels = combined)))

  ## ============ Compute minimal junction read coverage ===========
  if(verbose) message("Step 7: Computing minimal junction read coverage")

  ##  Add columns with the number of reads supporting the left and right splice
  ##  junction and the minimum of both
  key_sj <- with(sj, paste(seqnames, start - 1, end + 1, strand, sep = "."))
  key_ne_l <- paste(novel_exons$seqnames, novel_exons$lend, novel_exons$start,
                    novel_exons$strand, sep = ".")
  key_ne_r <- paste(novel_exons$seqnames, novel_exons$end, novel_exons$rstart,
                    novel_exons$strand, sep = ".")

  novel_exons$unique_left <- sj$unique[match(key_ne_l, key_sj)]
  novel_exons$unique_right <- sj$unique[match(key_ne_r, key_sj)]

  ## if the exon is a cassette exons, but one of the junctions is not in
  #SJ.out.tab (count = NA), # remove the prediction, because it is most likely
  #wrong (this only happens if we predict exons from reads) novel_exons <-
  #novel_exons %>% filter(!(!is.na(lend) & is.na(unique_left))) %>%
  #filter(!(!is.na(rstart) & is.na(unique_right)))

  ## take the minimum read coverage of both junctions
  novel_exons$min_reads <- pmin(novel_exons$unique_left,
                                novel_exons$unique_right, na.rm = TRUE)
  novel_exons
}


