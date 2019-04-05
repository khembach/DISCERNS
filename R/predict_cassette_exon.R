
#' Predict novel cassette exons
#'
#' Novel cassette exons are predicted from pairs of novel splice junctions (SJ)
#' that are located within an annotated intron and share the start and end
#' coordinates of the intron:
#' \preformatted{
#' X---------X   annotated intron
#' x---x         novel splice junction
#'       x---x   novel splice junction
#' X---NNN---X  predicted cassette exon (N)
#' }
#' First, the novel SJs are filtered: Only SJs that are located within an annotated
#' intron and that share their start or end coordinates with the intron are
#' retained. All introns that share both their start and end with a novel SJ are
#' tested for cassette exons. If the two novel SJs within an intron do not
#' overlap and are on the same strand, a novel cassette exon is predicted.
#'
#' @param sj_unann GRanges object with unannotated splice junctions.
#' @param introns GRanges object with intron annotations, e.g. the "introns" slot
#'  from the [prepare_annotation()] return object.
#'
#' @return List with two slots: "ne" is a data.frame with the coordinates of
#' the identifies cassette exons. It has 6 columns: `seqnames`, `lend`, `start`,
#'  `end`, `rstart` and `strand`.. "sj" is a GRanges object with the SJs from 
#' sj_unann that were not used to predict a cassette exons.
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges subsetByOverlaps
#' @export
#'
predict_cassette_exon <- function(sj_unann, introns) {
  within_ind <- unique(queryHits(findOverlaps(sj_unann, introns,
                                              type = "within")))
  within <- sj_unann[within_ind, ]

  start_hit <- findOverlaps(within, introns, type = "start")
  end_hit <- findOverlaps(within, introns, type = "end")
  ## filter out introns with only start/end
  start_hit <- start_hit[subjectHits(start_hit) %in% subjectHits(end_hit), ]
  end_hit <- end_hit[subjectHits(end_hit) %in% subjectHits(start_hit), ]

  ## filter out introns with more than one start or end junction
  start_intr_id <- table(subjectHits(start_hit))
  start_intr_id <- names(start_intr_id[start_intr_id == 1])
  end_intr_id <- table(subjectHits(end_hit))
  end_intr_id <- names(end_intr_id[end_intr_id == 1])

  intr_id<- as.integer(start_intr_id[start_intr_id %in% end_intr_id])

  start_hit <- start_hit[subjectHits(start_hit) %in% intr_id]
  starts <- within[queryHits(start_hit)[order(subjectHits(start_hit))]]

  end_hit <-  end_hit[subjectHits(end_hit) %in% intr_id]
  ends <- within[queryHits(end_hit)[order(subjectHits(end_hit))]]

  ## filter out overlapping SJs or SJs on different strands
  paired_junc_id <- which(end(starts) < start(ends) &
                            strand(starts) == strand(ends))
  starts <- starts[paired_junc_id]
  ends <- ends[paired_junc_id]

  if (length(starts) > 0 & length(ends) > 0) {
    novel_exons <- data.frame(seqnames = as.vector(seqnames(starts)),
                              lend = start(starts) - 1L,
                              start = end(starts) + 1L, end = start(ends) - 1L,
                              rstart = end(ends) + 1L,
                              strand = as.vector(strand(starts)),
                              stringsAsFactors = TRUE )
    sj_unann <- subsetByOverlaps(sj_unann, c(starts, ends),
                                 type = "equal", invert = TRUE)
  } else {
    novel_exons <- data.frame(seqnames = character(), lend = integer(),
                              start = integer(), end = integer(),
                              rstart = integer(), strand = character(),
                              stringsAsFactors = TRUE)
  }
  return(list(ne = novel_exons, sj = sj_unann))
}
