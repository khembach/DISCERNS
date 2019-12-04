#' Filter novel splice junctions based on overlapping genes
#'
#' This function filters splice junctions based on the genes that overlap with
#' the start and end of the splice junction, and that are crossed by the splice
#' junction.
#'
#' @param sg Integer vector. Index of the genes that overlap with the SJ start.
#' @param eg Integer vector. Index of the genes that overlap with the SJ end.
#' @param sjg Integer vector. Index of the genes that overlap with the SJ.
#'
#' @return Logical scalar. FALSE if the start and end of the SJ overlap with
#'   different genes.
#' @export
#'
filter_trans_splice_junctions <- function(sg, eg, sjg) {
  if (!is.null(sg)) { ## sg is defined
    if (!is.null(eg)){  ## eg is defined
      if (setequal(sg, eg)) { ## both are defined and equal
        TRUE          
      } else {
        FALSE 
      }
    } else {  ## sg is defined and eg is NULL
      if (setequal(c(sg, eg), sjg)) {  ## start and end cover the same genes
        TRUE
      } else {
        FALSE   
      }
    }
  } else if(!is.null(eg)) { ## sg is not defined and eg is defined
    if (setequal(c(sg, eg), sjg)) {  ## start and end cover the same genes
      TRUE
    } else {
      FALSE   
    }
  } else {  ## both are not defined
    FALSE
  }
}


#' Determine which side of the splice junction touches an exon
#'
#' Find all exons with the same strand as the splice junction and that touch it.
#' If a sj touches two exons, the exons must come from the same gene.
#' @param exons GRanges object. Annotated exons.
#' @param seqnames Character string or factor. Seqnames (chromosome) of the
#'   splice junction.
#' @param start Integer scalar. Start of the splice junction.
#' @param end Integer scalar. End of the splice junction.
#' @param strand Character string or factor. Strand of the splice junction ("-"
#'   or "+").
#'
#' @return A string that defines which side of the splice junction is touching
#'   an exon: "both", "start", "end", NA if the splice junction does not touch
#'   an exon or if the exons come from different genes.
#' @export
#'
sj_touching_exon <- function(seqnames, start, end, strand, exons) {
  ex <- exons[seqnames(exons) == seqnames]
  
  s <- which(start - 1 == end(ex))
  s <- ex[s]   ## all touching exons at start of sj
  s <- s[strand(s) == strand, ] ## all touching exons with same strand
  
  e <- which((end + 1) == start(ex))
  e <- ex[e]   ## all touching exons at end of sj
  e <- e[strand(e) == strand, ]
  
  if (length(s) > 0){
    if (length(e) > 0){
      ## pair of touching exons from the same gene
      if (any(mcols(s)$gene_id %in% mcols(e)$gene_id)) {
        "both"
      } else NA
    } else "start"
  } else if (length(e) > 0){
    "end"
  } else NA
}


#' Match novel splice junctions within an intron
#'
#' @param s GRanges object. Junctions from the start of the intron.
#' @param e GRanges object. Junctions from the end of the intron.
#'
#' @return A vector with the seqnames, end of the preceding exon, start of the
#'   novel exon, end of the novel exon, start of the consecutive exon, and the
#'   strand.
#' @export
#'
match_sj_in_intron <- function(s, e) {
  if (all( length(s) == 1, length(e) == 1,
           isDisjoint(c(s, e), ignore.strand=FALSE))) {
    return(c(as.vector(seqnames(s)), start(s) - 1, end(s) + 1,
             start(e) - 1, end(e) + 1, as.vector(strand(s))))
  } ## TODO: more than two novel splice junctions per intron
}


#' Transcript range
#'
#' This function returns a GRanges object with the start and end of a transcript
#' and takes a GRanges object with the transcript's exons as input.
#'
#' @param gr GRanges object. All exons of a transcript.
#'
#' @return GRanges object with the start and end of the transcript.
#' @export
#' 
transcript_range <- function(gr) {
  ## TODO: This is only needed for the simulated data, because the transcript
  ## annotation still contains the removed exons in the GTF file. For real data,
  ## we can simply take to coordinates of the entry with type="transcript".

  GRanges(seqnames = seqnames(gr)[1],
          ranges = IRanges(min(start(gr)), max(end(gr))),
          strand = strand(gr)[1])
}


#' Determine terminal end of a splice junction
#'
#' This function computes if a splice junction (SJ) splices to a terminal exons.
#' If the exon at the start of the SJ is terminal, the function returns "start",
#' "end" if it is the exon at the end and "NA" if none of the exons are
#' terminal.
#'
#' @param j GRanges object. Novel splice junction.
#' @inheritParams get_second_sj
#'
#' @return Character string; either "start", "end" or "NA".
#'
#' @importFrom GenomicFeatures transcripts
#' @export
#'
#' @section Details:
#'   The function first overlaps the SJ with all genes and then with all
#'   transcripts of the overlapping gene.
#'
which_exon_terminal <- function(j, txdb = txdb, gtxdb = gtxdb, ebyTr = ebyTr) {
  ### TODO: maybe we have to include the position +-1 of the sj, in case of
  ### terminal sj that is outside of the annotated gene boundary ?

  ## all genes that overlap with the SJ
  genes <- mcols(subsetByOverlaps(gtxdb, j))$gene_id
  tr <- transcripts(txdb, filter = list(gene_id = genes))
  tr <- mcols(tr)$tx_name  ## all transcripts from the genes
  e <- ebyTr[tr]
  ## the range of all transcripts
  tr_range <- unlist(GRangesList(lapply(e, transcript_range)))
  ## all transcripts that do NOT overlap with the sj
  tr <- names(subsetByOverlaps(tr_range, j, invert = TRUE))

  ## TODO: Use this for real data, where the transcript boundaries correspond to
  #the start and end of the first and last exon, but it does not work for the
  #simulated data, because the "transcript" entry in the gtf file still contains
  #the removed exons!
  #tr <- subsetByOverlaps(transcripts(txdb, filter = list(gene_id = genes)),
  #                                   j, invert = TRUE)

  e <- unlist(ebyTr[tr]) ## all exons from the non-overlapping transcripts
  if (length(e) == 0) {  # Sj overlaps with all transcripts of the gene.
    return(NA)
  }
  ## does the sj touch any of the exons?
  start <- start(j) - 1 == end(e)
  end <- end(j) + 1 == start(e)
  return(ifelse(any(start), "start", ifelse(any(end), "end", "NA")))
}
