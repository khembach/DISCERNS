#' Junction reads from BAM file
#'
#' The function takes a BAM file as input and returns all reads with "N" in the
#' CIGAR string and that are properly paired.
#'
#' The `yield_size` param determines the runtime: The bigger, the faster. If
#' possible, use at least 200000. The same goes for the `tile_width` param: the
#' bigger the faster. Also, the BAM file can be read in parallel (with mclapply)
#' if `cores` > 1.
#' 
#' Per tile, we only keep the reads that have their end location inside the
#' tile. This prevents having duplicate reads in the output, in case the read
#' overlaps the tile boundary.
#' 
#' Note: The human gneome has ~3 billion base pairs. If we have a BAM file with
#' 100 million reads --> assuming uniform read coverage: 0.033 reads per bp
#' a genome tile of 1e7 contains 333k reads.
#'
#' @inheritParams find_novel_exons
#'
#'@return GAlignments object with all junction reads from the `bam` file.
#'@importFrom Rsamtools BamFile scanBamFlag bamFlagAsBitMatrix 
#'@importFrom GenomicAlignments readGAlignments cigarOpTable cigar
#'@importFrom GenomicFiles reduceByYield
#'@importFrom parallel mclapply
#'
#'@export
#'
filter_junction_reads <- function(bam, lib_type = "PE", 
                                       stranded = "reverse", cores = 1, 
                                       yield_size = 200000, tile_width = 1e7) {
  if (!lib_type %in% c("SE", "PE")) {
    stop('Parameter lib_type has to be either "SE" or "PE".')
  }
  if (!stranded %in% c("unstranded", "forward", "reverse")) {
    stop('Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  }
  
  bf <- BamFile(bam, yieldSize = yield_size)
  tiles <- tileGenome(seqinfo(bf), tilewidth = tile_width)
  
  ## Read and filter function definitions
  read_bam_PE <- function(bf, tile, ...) {
    flag0 <- scanBamFlag(isProperPair = TRUE)
    param <- ScanBamParam(what = c("qname", "flag"), flag = flag0, which = tile)
    reads <- readGAlignments(bf, param = param)
    ## make sure that all reads are within the tile, sometimes pairs might be
    ## located on different tiles
    reads[end(reads) <= end(tile)]
  }
  
  read_bam_SE <- function(bf, tile, ...) {
    param <- ScanBamParam(what = c("qname", "flag"), which = tile)
    reads <- readGAlignments(bf, param = param)
    reads[end(reads) <= end(tile)]
  }
  
  filter_bam_PE <- function(reads, stranded, ...) {
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
  
  filter_bam_SE <- function(reads, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag, bitnames = c("isMinusStrand"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    strand(reads) <- ifelse(mcols(reads)$isMinusStrand == 0, "-", "+")
    reads
  }
  
  filter_bam_UNSTRANDED <- function(reads, ...) {
    reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
  }
 
  if (lib_type == "PE") {
    if( stranded == "unstranded") {
      mclapply(tiles, function(tile) {
        reads <- read_bam_PE(bf, tile)
        filter_bam_UNSTRANDED(reads) }, mc.cores = cores)
    } else {
      mclapply(tiles, function(tile) {
        reads <- read_bam_PE(bf, tile)
        filter_bam_PE(reads, stranded) }, mc.cores = cores)
    }
  } else if (stranded != "reverse") { # SE
    mclapply(tiles, function(tile) {
      reads <- read_bam_SE(bf, tile)
      filter_bam_UNSTRANDED(reads) }, mc.cores = cores)
  } else {
    mclapply(tiles, function(tile) {
      reads <- read_bam_SE(bf, tile)
      filter_bam_SE(reads) }, mc.cores = cores)
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
#' @param junc_reads GAlignments object with junction reads.
#' @inheritParams find_novel_exons
#' 
#' @return data.table with the coordinates of the predicted novel exon. It has 6
#'   columns: seqnames, lend, start, end, rstart and strand
#'
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom data.table data.table as.data.table setkeyv setorder %chin%
#'   setnames
#'
#' @export
#'
predict_jr_exon <- function(junc_reads, annotation) {
  if (length(junc_reads) == 0){
    stop("The GAlignments object for parameter junc_reads is empty.")
  }
  
  ## Keep all reads with 2 junctions and also look for novel junction
  ## combinations. Filter all reads with 2 "N" in cigar
  two_jr <- junc_reads[njunc(junc_reads) == 2, ]
  cigar_junc <- cigarRangesAlongReferenceSpace(cigar = cigar(two_jr),
                                               flag = mcols(two_jr)$flag,
                                               pos = start(two_jr),
                                               ops = "N")
  names(cigar_junc) <- 1:length(two_jr)
  cigar_junc <- unlist(cigar_junc)
  cigar_junc <- GRanges(seqnames(two_jr)[as.integer(names(cigar_junc))],
                           cigar_junc,
                           strand(two_jr)[as.integer(names(cigar_junc))])
  rm(two_jr)
  ## we join the read junctions with the annotated introns: start is the first
  ## nucleotide in the intron and end the last nucleotide (excluding exons)
  intr_gtf <- unlist(intronsByTranscript(annotation[["txdb"]], 
                                         use.names = TRUE))
  ## TODO: do we need this? Are there transcripts with duplicate introns?
  # intr_gtf$transcript_id <- names(intr_gtf)
  # intr_idx <- paste0(start(intr_gtf), ":", end(intr_gtf), "_",
  #                    mcols(intr_gtf)$transcript_id)
  # intr_gtf<- intr_gtf[match(unique(intr_idx), intr_idx)]
  ##only keep the unique junction-transcript_id combinations
  
  ## we compare the read junctions to annotated transcript introns
  olap <- findOverlaps(cigar_junc, intr_gtf, type = "equal")
  
  tmp <- data.table(read = names(cigar_junc)[queryHits(olap)],
                    tr = names(intr_gtf[subjectHits(olap)]))
  setkeyv(tmp, c("read", "tr"))
  ## all reads where both junctions are annotated in a transcript
  read_id <- tmp[duplicated(tmp, by = c("read", "tr")), read]
  rm(tmp, olap)

  ## all reads with junctions that are not annotated in a transcripts
  cigar_junc_pred <- cigar_junc[!names(cigar_junc) %chin% read_id, ]
  
  ## All novel reads, ordered by start position
  novel_reads <- as.data.table(cigar_junc_pred) 
  novel_reads[, names := names(cigar_junc_pred)]
  setorder(novel_reads, "start")

  ## Per read, compute the coordinates of the novel exon
  read_pred <- novel_reads[, by = names,
                           .(lend = start[1] - 1L, start1 = end[1] + 1L, 
                             end1 = start[2] - 1L, rstart = end[2] + 1L, 
                             seqnames = unique(seqnames), 
                             strand = unique(strand))]
 
  read_pred <- read_pred[,!"names"]
  setnames(read_pred, old = c("start1", "end1"), new = c("start", "end"))
  read_pred <- unique(read_pred)
  
  ## keep all predictions where the exon itself is not annotated
  read_pred <- as.data.table(subsetByOverlaps(
    GRanges(read_pred), annotation[["exons"]], type = "equal", invert = TRUE))
  read_pred <- read_pred[,!"width"]
  
  ### Keep all exons that are located within gene boundaries
  read_pred <- read_pred[
    unique(queryHits(findOverlaps(
      GRanges(read_pred), genes(annotation[["txdb"]]), type = "within"))), ]
  read_pred
}



#' Predict novel exons from read pairs with two splice junctions
#'
#' Novel exons are predicted from paired-end reads where each read spans one
#' splice junction. First, the read pairs are filtered and the distance between
#' the end of the first junction and the start of the second junction has to be
#' `< 2 * (read_length - overhang_min) + min_intron_size`. Here, `read_length`
#' is the length of the reads, `overhang_min` is the minmal required read
#' overhang over a splice junction of the alignment tool and `min_intron_size`
#' is the minimal required intron length of the alignment tool. For example,
#' paired-end reads with a lenght of 101 nts and a minimal overhang of 6 and a
#' minimal intron length of 21 allow a distance of at most 211 nucleotides
#' between the two splice junctions: `2 * (101 - 6) + 21 = 211`. If the distance
#' between the two splice junctions exceeds the limit, it cannot be guaranteed
#' that the junctions are connected to the same exon. Splice junction pairs that
#' are already annotated in a transcript are removed. Novel exons are predicted
#' from the remaining splice junction pairs.
#' 
#' @param junc_reads GAlignments object with junction reads.
#' @inheritParams find_novel_exons
#'
#' @return data.frame with the coordinates of the predicted novel exon. It has 6
#'   columns: `seqnames`, `lend`, `start`, `end`, `rstart` and `strand`.
#'
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom data.table as.data.table := setkeyv %chin%
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse
#'
#' @export
#'
predict_jrp_exon <- function(junc_reads, annotation, 
                             read_length = 101, overhang_min = 12, 
                             min_intron_size = 21) {
  if (any(read_length <= 0, overhang_min <= 0, min_intron_size <= 0)) {
    stop('Parameters "read_length", "overhang_min" and "min_intron_size" must be >= 0.')
  }
  
  if (length(junc_reads) == 0){
    stop("The GAlignments object for parameter junc_reads is empty.")
  }
  
  junc_rp <- junc_reads[njunc(junc_reads) == 1]
 
  ## keep all reads with duplicate read ids (both reads in pair have a sj)
  junc_rp <- junc_rp[mcols(junc_rp)$qname %chin% 
                       mcols(junc_rp)$qname[duplicated(mcols(junc_rp)$qname)]]
  
  ## keep all reads pairs where the junctions are already annotated
  cigar_jp <- cigarRangesAlongReferenceSpace(cigar = cigar(junc_rp), 
                                             flag = mcols(junc_rp)$flag,
                                             pos = start(junc_rp),
                                             ops = "N")
  ## add strand and seqnames info 
  names(cigar_jp) <- mcols(junc_rp)$qname
  cigar_jp <- as.data.table(unlist(cigar_jp))
  cigar_jp <- cigar_jp[, !"width"]
  m <- match(cigar_jp$names, mcols(junc_rp)$qname)
  cigar_jp[ , ':=' (seqnames = as.character(seqnames(junc_rp)[m]),
                    strand = as.character(strand(junc_rp)[m]))]
  rm(junc_rp, m)

  ## remove all duplicate junctions within one read pair and only keep the read
  ## pairs with two different junctions
  
  ## get unique reads, count read pairs with 2 splice junctions
  ## sort the table by read names
  setkeyv(cigar_jp, "names")  ## sorting

  ## Reads might be duplicated because they overlap the tile boundaries. A read
  ## name can appear at most 4 times - 4 consecutive rows - if both reads are
  ## duplicated.
  ## This function consecutively processes all reads. We remove all duplicated
  ## reads and all read pairs with only one distince splice junction (both reads
  ## in the pair have the same SJ).
  ind <- get_unique_ind(cigar_jp)
  cigar_jp <- cigar_jp[ind]

  ## compute length of the novel exon and discard all exons that are longer than
  ## the minimal intron length in STAR
  max_pairs_exon_len <- 2*(read_length - overhang_min) + min_intron_size
  
  ## This function computes the lenght of the novel exon and removes the
  ## prediction if the novel exon is longer than `max_pairs_exon_len`
  
  x <- filter_exon_length(cigar_jp, max_pairs_exon_len)
  x <- cigar_jp[x, names]
  
  cigar_jp <- cigar_jp[cigar_jp$names %chin% x, ]
  cigar_jp_gr <- GRanges(cigar_jp)
  
  inbytx <- intronsByTranscript(annotation[["txdb"]], use.names = TRUE)
  intr_gtf <- unlist(inbytx)
  
  ## unify the seqnames of the novel SJs and the annotation
  intr_gtf <- intr_gtf[seqnames(intr_gtf) %in% seqlevels(cigar_jp_gr)]
  GenomeInfoDb::seqlevels(intr_gtf) <- seqlevelsInUse(intr_gtf)
  
  ## Keep all reads with annotated junctions
  ## names of all reads with unannotated junctions
  x <- mcols(cigar_jp_gr)$names[-unique(
    queryHits(findOverlaps(cigar_jp_gr,intr_gtf,type = "equal")))]
  ## all read pairs with only annotated junctions
  cigar_jp_gr <- cigar_jp_gr[!mcols(cigar_jp_gr)$names %chin% x ]
  ## we compare the read junctions to annotated transcript introns
  olap <- findOverlaps(cigar_jp_gr, intr_gtf, type = "equal")
  
  tmp <- data.table(read = mcols(cigar_jp_gr)$names[queryHits(olap)],
                    tr = names(intr_gtf[subjectHits(olap)]))
  setkeyv(tmp, c("read", "tr"))
  
  ## all reads where both junctions are annotated in a transcript
  read_id <- tmp[duplicated(tmp, by = c("read", "tr")), read]
  rm(tmp, olap)
  
  novel_reads <- as.data.table(
    cigar_jp_gr[mcols(cigar_jp_gr)$names %chin% read_id, ])
  setkeyv(novel_reads, c("names", "start"))
  
  
  ## This function computes the coordinates of the novel exons. The input
  ## data.table is sorted by read name and SJ start position.
  
  ## Predict the coordinates of the novel exons
  read_pairs_pred <- pred_exon_coord(novel_reads)
  
  ## Keep all predictions, where the exon itself is not jet annotated
  read_pairs_pred <- as.data.table(subsetByOverlaps(GRanges(read_pairs_pred),
                                                    annotation[["exons"]],
                                                    type = "equal",
                                                    invert = TRUE))
  read_pairs_pred <- read_pairs_pred[,!"width"]
  
  
  ### Keep all exons that are located within gene boundaries (both annotated
  ### junctions are from the same gene)
  read_pairs_pred <- read_pairs_pred[unique(queryHits(findOverlaps(
    GRanges(read_pairs_pred[, seqnames],
            IRanges(read_pairs_pred[, lend], read_pairs_pred[, rstart]),
            read_pairs_pred[,strand]),
    genes(annotation[["txdb"]]), type = "within"))), ]
  read_pairs_pred
}
