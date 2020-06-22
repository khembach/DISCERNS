#' Prepare anntation
#'
#' Read exon and intron annotation from a GTF file.
#'
#' @param gtf_file Character string. Path to a GTF file.
#'
#' @return List with exon and intron annotations as GRanges and the `TxDB`
#'   object. The slots are `exons`, `introns` and `txdb`.
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicFeatures makeTxDbFromGRanges intronsByTranscript
#'
#' @export
#' 
#' @examples
#' gtf <- system.file("extdata", "selected.gtf", package = "DISCERNS", 
#'                    mustWork = TRUE)
#' anno <- prepare_annotation(gtf)
#' names(anno)
prepare_annotation <- function(gtf_file) {
  if(!file.exists(gtf_file)){
    stop( paste0("File ", gtf_file, " does not exist."))
  }
  
  gtf <- import(gtf_file)
  exons <- gtf[mcols(gtf)$type =="exon", ]
  
  txdb <- makeTxDbFromGRanges(exons)
  inbytx <- intronsByTranscript(txdb, use.names=TRUE)
  introns <- unique(unlist(inbytx))
  
  list(exons = exons, introns = introns, txdb = txdb)
}
