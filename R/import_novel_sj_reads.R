

#' Import reads with novel splice junctions
#'
#' This functions takes a BAM file and a GRanges object with unannotated splice
#' junctions as input and returns a GAlignments object with all reads from the
#' BAM file that contain any of the novel splice junctions.
#'
#' @param sj_unann GRanges object with unannotated splice junctions.
#' @inheritParams find_novel_exons
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



#' Import splice junction reads in a given genomic region
#'
#' This function imports all splice junction reads within a specified genomic
#' region.
#'
#' @param bf BamFile object.
#' @param g GRanges object specifying the genomic region from which to import
#'   the reads.
#' @param region
#'
#' @return GAlignments object with all reads from `bf` overlapping `i`.
#' @export
#' 
get_reads <- function(bf, g, region){
  reads <- readGAlignments(bf, index = bam, with.which_label = TRUE,
                           param = ScanBamParam(which = g, what = c("qname")))
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
  
  reads[df$read_nr[df$junction == df$true_junction]]
}

#' Import reads with novel splice junctions in parallel
#'
#' This functions takes a BAM file and a GRanges object with unannotated splice
#' junctions as input and returns a GAlignments object with all reads from the
#' BAM file that contain any of the novel splice junctions. The BAM files is
#' read in parallel with mclapply using the specified number of `cores`.
#'
#' @param sj_unann GRanges object with unannotated splice junctions.
#' @inheritParams find_novel_exons
#'
#' @return GAlignments object with all reads from the BAM file that contain a
#'   novel splice junctions.
#'
#' @importFrom Rsamtools ScanBamParam BamFile
#' @importFrom GenomicAlignments readGAlignments njunc junctions
#' @importFrom parallel mclapply
#'
#' @export
#' 
import_novel_sj_reads_parallel <- function(bam, yield_size = 200000, sj_unann, 
                                           cores = 1) {
  bf <- BamFile(bam, yieldSize = yield_size)
   mclapply(seq_along(sj_unann), function(i) 
    get_reads(bf, sj_unann[i]), mc.cores = cores)
}

