#' Prepare exon and intron annotation from GTF file
#'
#' @param GTF Path to GTF file
#'
#' @return List with exon and intron annotations as GRanges and the TxDB object.
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicFeatures makeTxDbFromGRanges intronsByTranscript
#'
#' @export
#'
prepare_annotation <- function(GTF) {

  gtf <- import(GTF)
  exons <- gtf[mcols(gtf)$type =="exon", ]

  txdb <- makeTxDbFromGRanges(exons)
  inbytx <- intronsByTranscript(txdb, use.names=TRUE)
  introns <- unique(unlist(inbytx))

  list(exons = exons, introns = introns, txdb = txdb)
}
