

#' Import reads with novel splice junctions
#'
#' This functions takes a BAM file and a GRanges object with unannotated splice
#' junctions as input and returns a GAlignments object with all reads from the
#' BAM file that contain any of the novel splice junctions.
#'
#' @param bam The path to the BAM file.
#' @param sj_unann GRanges object with unannotated splice junctions.
#'
#' @return GAlignments object with all reads from the BAM file that contain a
#'   novel splice junctions.
#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments njunc junctions
#'
#' @export
#'
import_novel_sj_reads <- function(bam, sj_unann) {

  param <- ScanBamParam(which = sj_unann, what = c("qname"))
  reads <- readGAlignments(bam, index = bam, with.which_label = TRUE,
                           param = param)
  reads <- reads[njunc(reads) > 0, ]
  junc <- junctions(reads, use.mcols = TRUE)
  names(junc) <- 1:length(junc)
  true_junc <- as.character(mcols(junc)$which_label)
  junc <- unlist(junc)
  df <- data.frame(read_nr = as.integer(names(junc)),
                   junction = paste0(seqnames(junc), ":", start(junc),
                                     "-", end(junc)),
                   true_junction = rep(true_junc, njunc(reads)),
                   stringsAsFactors = FALSE)

  reads <- reads[df$read_nr[df$junction == df$true_junction]]
  reads
}
