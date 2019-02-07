#' Filter reads with splice junctions from BAM file
#'
#' The function Filters all reads with "N" in the CIGAR string and that are
#' properly paired (-f 2) from a BAM file.
#'
#' @param bam The path to the BAM file.
#'
#' @return GAlignments object with junction reads
#'
#' @importFrom data.table fread
#' @importFrom Rsamtools bamFlagAsBitMatrix
#' @export
#'
filter_junction_reads <- function(bam){
  junc_reads <- fread(cmd = paste0("samtools view  -f 2 ", BAM,
                           " | awk '{if ($6 ~ /N/) print $1, $2, $3, $4, $6}'"))
  names(junc_reads) <- c("qname", "flag", "seqnames", "pos", "cigar")

  ## We infer the read strand from the flag: our data comes from Illumina HiSeq
  # 2000 (stranded TruSeq libary preparation with dUTPs ) Thus, the the last
  # read in pair determines the strand of the junction --> reverse the strand
  # of the first reads
  junc_reads <- cbind(junc_reads,
                      bamFlagAsBitMatrix(junc_reads$flag,
                                         bitnames=c("isMinusStrand",
                                                    "isFirstMateRead")) )

  ## TODO: make sure that the strand of both reads is correct --> use a
  ## parameter for the type of sequencing

  ## isMinusStrand==0 & isFirstMateRead==0 --> "+"
  ## isMinusStrand==0 & isFirstMateRead==1 --> "-"
  ## isMinusStrand==1 & isFirstMateRead==0 --> "-"
  ## isMinusStrand==1 & isFirstMateRead==1 --> "+"
  junc_reads[, strand := ifelse(isMinusStrand + isFirstMateRead == 0, "+",
                                ifelse(isMinusStrand + isFirstMateRead == 1, "-",
                                       "+"))]
  junc_reads
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
#' @param junc_reads GAlignments object with junction read
#' @param txdb TxDb object with gene annotations
#'
#' @return data.frame with the coordinates of the predicted novel exon. It has 6
#'   columns: seqnames, lend, start, end, rstart and strand
#'
#' @importFrom stringr str_count
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom dplyr group_by
#'
#' @export
#'
predict_jr_exon <- function(junc_reads, annotation){
  ##  keep all reads with 2 junctions and also look for novel junction combinations
  ## filter all reads with 2 "N" in cigar
  two_junc_reads<- junc_reads[str_count(junc_reads$cigar, "N") == 2, ]
  two_junc_reads <- two_junc_reads[!duplicated(two_junc_reads[, -"qname",
                                                              with = FALSE]), ]

  cigar_junc <- cigarRangesAlongReferenceSpace(cigar = two_junc_reads$cigar,
                                               flag = two_junc_reads$flag,
                                               pos = two_junc_reads$pos,
                                               ops = "N")
  names(cigar_junc) <- 1:nrow(two_junc_reads)
  cigar_junc_gr <- unlist(cigar_junc)
  cigar_junc_gr <- GRanges(two_junc_reads$seqnames[as.integer(
                                                        names(cigar_junc_gr))],
                           cigar_junc_gr,
                           two_junc_reads$strand[as.integer(
                                                        names(cigar_junc_gr))])

  ## we join the read junctions with the annotated introns start is the first
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
  olap <- findOverlaps(cigar_junc_gr, intr_gtf, type ="equal")

  ## test is TRUE if there is a transcript that contains both junctions
  tmp <- data.frame(read = names(cigar_junc_gr)[queryHits(olap)] ,
                    tr = names(intr_gtf[subjectHits(olap)])) %>%
    group_by(read,tr) %>%
    count(.) %>%
    group_by(read) %>%
    summarise(test = any(n == 2))

  read_id <- tmp %>%
    filter(test==T) %>%
    pull(read) ## all reads that are already annotated

  ## all reads with junctions that are not annotated in a transcripts
  cigar_junc_gr_pred <- cigar_junc_gr[!names(cigar_junc_gr) %in% read_id, ]

  novel_reads <-  as.data.table(cigar_junc_gr_pred) %>%
    mutate(names = names(cigar_junc_gr_pred))
  novel_reads <- novel_reads[order(novel_reads$start), ]
  read_pred <- novel_reads %>%
    group_by(names) %>%
    dplyr::summarise(lend = nth(start-1, 1),
                     start1 = nth(end+1, 1),
                     end1 = nth(start-1, 2),
                     rstart = nth(end+1, 2),
                     seqnames = unique(seqnames), strand = unique(strand))

  ## TODO: only keep predictions with >- x supporting reads
  read_pred <- read_pred %>%
    select(-names) %>%
    rename(start = start1, end = end1) %>%
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
