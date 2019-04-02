#' Junction reads from BAM file
#'
#' The function takes a BAM file as input and returns all reads with "N" in the
#' CIGAR string and that are properly paired.
#'
#' The yieldSize param determines the runtime: The bigger, the faster. If possible, use at
#' least 200000.
#'
#'@param bam The path to the BAM file.
#'@param yield_size Integer scalar. The number of reads that should be read in
#'  each chunk.
#'@param lib_type Character scalar. Type of the sequencing library: either "SE"
#'  (single-end) or "PE" (paired-end). Default: "PE"
#'@param stranded Character scalar. Strand type of the sequencing protocol:
#'  "unstranded" for unstranded protocols; "forward" or "reverse" for stranded
#'  protocols. "forward" refers to protocols where the first read comes from the
#'  forward (sense) strand. Default: "reverse"
#'
#'@return GAlignments object with all junction reads from the bam file.
#'@importFrom Rsamtools BamFile scanBamFlag bamFlagAsBitMatrix 
#'@importFrom GenomicAlignments readGAlignments cigarOpTable cigar
#'@importFrom GenomicFiles reduceByYield
#'
#'@export
#'
filter_junction_reads <- function(bam, yield_size = 200000, 
                                  lib_type = "PE", stranded = "reverse") {
  if (!lib_type %in% c("SE", "PE")) {
    stop('Parameter lib_type has to be either "SE" or "PE".')
  }
  if (!stranded %in% c("unstranded", "forward", "reverse")) {
    stop('Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  }
  bf <- BamFile(bam, yieldSize = yield_size)
  
  YIELD <- function(x, ...) {
    flag0 <- scanBamFlag(isProperPair = TRUE)
    param <- ScanBamParam(what = c("qname", "flag"), flag = flag0)
    readGAlignments(x, param = param)
  }
  
  YIELD_SE <- function(x, ...) {
    param <- ScanBamParam(what = c("qname", "flag"))
    readGAlignments(x, param = param)
  }
  
  MAP <- function(reads, stranded, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag,
                               bitnames = c("isMinusStrand", "isFirstMateRead"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    mcols(reads)$isFirstMateRead <- bfbm[ , "isFirstMateRead"]
    
    if (stranded == "forward") {
      ## isMinusStrand==0 & isFirstMateRead==0 --> "-"
      ## isMinusStrand==0 & isFirstMateRead==1 --> "+"
      ## isMinusStrand==1 & isFirstMateRead==0 --> "+"
      ## isMinusStrand==1 & isFirstMateRead==1 --> "-" 
      strand(reads) <- ifelse(mcols(reads)$isMinusStrand +
                                mcols(reads)$isFirstMateRead == 0, "-",
                              ifelse(mcols(reads)$isMinusStrand +
                                       mcols(reads)$isFirstMateRead == 1, "+",
                                     "-"))
      reads
    } else {  # reverse
      ## isMinusStrand==0 & isFirstMateRead==0 --> "+"
      ## isMinusStrand==0 & isFirstMateRead==1 --> "-"
      ## isMinusStrand==1 & isFirstMateRead==0 --> "-"
      ## isMinusStrand==1 & isFirstMateRead==1 --> "+"
      strand(reads) <- ifelse(mcols(reads)$isMinusStrand +
                                mcols(reads)$isFirstMateRead == 0, "+",
                              ifelse(mcols(reads)$isMinusStrand +
                                       mcols(reads)$isFirstMateRead == 1, "-",
                                     "+"))
      reads
    }
  } 
  
  MAP_SE_REVERSE <- function(reads, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag, bitnames = c("isMinusStrand"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    strand(reads) <- ifelse(mcols(reads)$isMinusStrand == 0, "-", "+")
    reads
  }
  
  MAP_UNSTRANDED <- function(reads, ...) {
    reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
  }
  
  DONE <- function(value) {
    length(value) == 0L
  }
  
  if (lib_type == "PE") {
    if( stranded == "unstranded") {
      reduceByYield(bf, YIELD, MAP_UNSTRANDED, REDUCE = c, DONE, parallel = FALSE)
    } else {
      reduceByYield(bf, YIELD, MAP, REDUCE = c, DONE, parallel = FALSE, 
                    stranded = stranded)
    }
  } else if (stranded != "reverse") { # SE
    reduceByYield(bf, YIELD_SE, MAP_UNSTRANDED, REDUCE = c, DONE, 
                  parallel = FALSE)
  } else { 
    reduceByYield(bf, YIELD_SE, MAP_SE_REVERSE, REDUCE = c, DONE, 
                  parallel = FALSE)
  }
}


#' Predict novel exons from reads with 2 junctions
#'
#' Novel exons are predicted from novel combinations of annotated splice
#' junctions. Reads spanning two splice junctions are filtered from the set of
#' all spliced reads. Each of the splice junction pairs are compared with the
#' annotated splice junctions per transcript and already annotated splice
#' junction pairs are removed. In the last step, the predicted exon that is
#' defined by the two splice junctions is returned if it is contained within the
#' boundaries of any annotated gene. This prevents the prediction of false
#' positive exons from wrongly mapped reads.
#'
#' @param annotation List with exon and intron annotation as GRanges. Created
#'   with prepare_annotation().
#' @param junc_reads GAlignments object with junction read.
#'
#' @return data.frame with the coordinates of the predicted novel exon. It has 6
#'   columns: seqnames, lend, start, end, rstart and strand
#'
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom dplyr group_by summarise filter pull nth select rename ungroup n
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table
#' @importFrom rlang .data
#'
#' @export
#'
predict_jr_exon <- function(junc_reads, annotation) {
  ##  keep all reads with 2 junctions and also look for novel junction combinations
  ## filter all reads with 2 "N" in cigar
  two_jr <- junc_reads[njunc(junc_reads) == 2, ]
  cigar_junc <- cigarRangesAlongReferenceSpace(cigar = cigar(two_jr),
                                               flag = mcols(two_jr)$flag,
                                               pos = start(two_jr),
                                               ops = "N")
  names(cigar_junc) <- 1:length(two_jr)
  cigar_junc_gr <- unlist(cigar_junc)
  cigar_junc_gr <- GRanges(seqnames(two_jr)[as.integer(names(cigar_junc_gr))],
                           cigar_junc_gr,
                           strand(two_jr)[as.integer(names(cigar_junc_gr))])
  
  ## we join the read junctions with the annotated introns: start is the first
  ## nucleotide in the intron and end the last nucleotide (excluding exons)
  inbytx <- intronsByTranscript(annotation[["txdb"]], use.names = TRUE)
  intr_gtf <- unlist(inbytx)
  ## TODO: do we need this? Are there transcripts with duplicate introns?
  # intr_gtf$transcript_id <- names(intr_gtf)
  # intr_idx <- paste0(start(intr_gtf), ":", end(intr_gtf), "_",
  #                    mcols(intr_gtf)$transcript_id)
  # intr_gtf<- intr_gtf[match(unique(intr_idx), intr_idx)]
  ##only keep the unique junction-transcript_id combinations
  
  ## we compare the read junctions to annotated transcript introns
  olap <- findOverlaps(cigar_junc_gr, intr_gtf, type = "equal")
  
  ## test is TRUE if there is a transcript that contains both junctions
  tmp <- data.frame(read = names(cigar_junc_gr)[queryHits(olap)],
                    tr = names(intr_gtf[subjectHits(olap)])) %>%
    group_by(.data$read, .data$tr) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(.data$read) %>%
    summarise(test = any(.data$n == 2))
  
  read_id <- tmp %>%
    filter(.data$test == TRUE) %>%
    pull(.data$read) ## all reads that are already annotated
  
  ## all reads with junctions that are not annotated in a transcripts
  cigar_junc_gr_pred <- cigar_junc_gr[!names(cigar_junc_gr) %in% read_id, ]
  
  novel_reads <-  as.data.table(cigar_junc_gr_pred) %>%
    mutate(names = names(cigar_junc_gr_pred))
  novel_reads <- novel_reads[order(novel_reads$start), ]
  read_pred <- novel_reads %>%
    group_by(names) %>%
    summarise(lend = nth(start - 1, 1), start1 = nth(end + 1, 1),
              end1 = nth(start - 1, 2), rstart = nth(end + 1, 2),
              seqnames = unique(seqnames), strand = unique(strand))
  
  ## TODO: only keep predictions with >- x supporting reads
  read_pred <- read_pred %>%
    select(-names) %>%
    rename(start = .data$start1, end = .data$end1) %>%
    unique()
  
  ## keep all predictions where the exon itself is not annotated
  read_pred <- as.data.frame(subsetByOverlaps(
    GRanges(read_pred), annotation[["exons"]], type = "equal", invert = TRUE))
  read_pred <- read_pred %>% select(-width)
  
  ### Keep all exons that are located within gene boundaries
  read_pred <- read_pred[
    unique(queryHits(findOverlaps(
      GRanges(read_pred$seqnames, IRanges(read_pred$lend, read_pred$rstart),
              read_pred$strand ),
      genes(annotation[["txdb"]]), type = "within"))), ]
  read_pred
}


#' Predict novel exons from read pairs with two splice junctions
#'
#' Novel exons are predicted from paired-end reads where each read spans one
#' splice junction. First, the read pairs are filtered and the distance between
#' the end of the first junction and the start of the second junction has to be
#' `< 2*(readlength-minOverhang) + minIntronSize`. Here, `readlength` is the length
#' of the reads, `minOverhang` is the minmal required read overhang over a splice
#' junction of the alignment tool and `minIntronSize` is the minimal required
#' intron length of the alignment tool. For example, paired-end reads with a
#' lenght of 101 nts and a minimal overhang of 6 and a minimal intron length of
#' 21 allow a distance of at most 211 nucleotides between the two splice
#' junctions: `2*(101-6) + 21 = 211`. If the distance between the two splice
#' junctions exceeds the limit, it cannot be guaranteed that the junctions are
#' connected to the same exon. Splice junction pairs that are already annotated
#' in a transcript are removed. Novel exons are predicted from the remaining
#' splice junction pairs.
#'
#' @param annotation List with exon and intron annotation as GRanges. Created
#'   with prepare_annotation().
#' @param junc_reads GAlignments object with junction read.
#'
#' @return data.frame with the coordinates of the predicted novel exon. It has 6
#'   columns: `seqnames`, `lend`, `start`, `end`, `rstart` and `strand`.
#'
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr distinct group_by summarise filter pull nth select rename
#'   ungroup n
#' @importFrom data.table as.data.table :=
#'
#' @export
#'
predict_jrp_exon <- function(junc_reads, annotation) {
  junc_rp <- junc_reads[njunc(junc_reads) == 1, ]
  junc_rp <- junc_rp[duplicated(mcols(junc_rp)$qname) | duplicated(mcols(junc_rp)$qname,
                                                                   fromLast = TRUE), ]
  
  ## keep all reads pairs where the junctions are already annotated
  cigar_jp <- cigarRangesAlongReferenceSpace(cigar = cigar(junc_rp), 
                                             flag = mcols(junc_rp)$flag,
                                             pos = start(junc_rp),
                                             ops = "N")
  
  names(cigar_jp) <- mcols(junc_rp)$qname
  cigar_jp <- as.data.table(unlist(cigar_jp))
  m <- match(cigar_jp$names, mcols(junc_rp)$qname)
  cigar_jp[ , ':=' (seqnames = as.character(seqnames(junc_rp)[m]),
                    strand = as.character(strand(junc_rp)[m]))]
  
  ## remove all duplicate junctions within one read pair and only keep the read
  ## pairs with two different junctions
  cigar_jp <- cigar_jp %>% distinct
  two_junc_pairs <- cigar_jp %>%
    group_by(.data$names) %>%
    summarise(n = n()) %>%
    filter(.data$n == 2) %>%
    pull(.data$names)
  cigar_jp <- cigar_jp %>% filter(.data$names %in% two_junc_pairs)
  
  ## TODO: make this flexible for other read lengths and parameters
  max_pairs_exon_len <- 211
  x <- cigar_jp %>%
    group_by(.data$names) %>%
    select(.data$names, .data$start, .data$end) %>%
    summarise(d = nth(.data$start - 1, 2) - nth(.data$end + 1, 1)) %>%
    filter(.data$d < max_pairs_exon_len) %>%
    pull(.data$names)
  cigar_jp <- cigar_jp[cigar_jp$names %in% x, ]
  cigar_jp_gr <- GRanges(cigar_jp)
  
  inbytx <- intronsByTranscript(annotation[["txdb"]], use.names = TRUE)
  intr_gtf <- unlist(inbytx)
  
  ## Keep all reads with annotated junctions
  ## names of all reads with unannotated junctions
  x <- mcols(cigar_jp_gr)$names[-unique(queryHits(findOverlaps(cigar_jp_gr,
                                                               intr_gtf,
                                                               type = "equal")))]
  ## all read pairs with only annotated junctions
  cigar_jp_gr <- cigar_jp_gr[!mcols(cigar_jp_gr)$names %in% x ]
  ## we compare the read junctions to annotated transcript introns
  olap <- findOverlaps(cigar_jp_gr, intr_gtf, type = "equal")
  ## test is TRUE if there is a transcript that contains both junctions
  tmp <- data.frame(read = mcols(cigar_jp_gr)$names[queryHits(olap)],
                    tr = names(intr_gtf[subjectHits(olap)])) %>%
    group_by(.data$read, .data$tr) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(.data$read) %>%
    summarise(test = any(.data$n == 2))
  
  read_id <- tmp %>%
    filter(.data$test == FALSE) %>%
    pull(.data$read)
  
  cigar_jp_gr_pred <- cigar_jp_gr[mcols(cigar_jp_gr)$names %in% read_id, ]
  novel_reads <- as.data.table(cigar_jp_gr_pred)
  
  ## predict the novel exons based on the novel junctions
  novel_reads <- novel_reads[order(novel_reads$start), ]
  read_pairs_pred <- novel_reads %>% group_by(names) %>%
    summarise(lend = nth(start - 1, 1), start1 = nth(end + 1, 1),
              end1 = nth(start - 1, 2), rstart = nth(end + 1, 2),
              seqnames = unique(seqnames), strand = unique(strand))
  read_pairs_pred <- read_pairs_pred %>%
    select(-names) %>%
    rename(start = .data$start1, end = .data$end1) %>%
    unique()
  
  ## Keep all predictions, where the exon itself is not jet annotated
  read_pairs_pred <- as.data.frame(subsetByOverlaps(GRanges(read_pairs_pred),
                                                    annotation[["exons"]],
                                                    type = "equal",
                                                    invert = TRUE))
  read_pairs_pred <- read_pairs_pred %>% select(-width)
  read_pairs_pred$seqnames <- read_pairs_pred$seqnames
  read_pairs_pred
}
