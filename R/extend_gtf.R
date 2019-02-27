#' Find all transcripts that contain a given exon
#'
#' Given a novel exon, find the IDs of all transcripts that contain the upstream
#' and downstream exons.
#'
#' @param exons GRanges object with exon annotations.
#' @param seqn Scalar, seqname of the novel exon.
#' @param lend Integer scalar, end of the upstream exon.
#' @param start Integer scalar, start of the novel exon.
#' @param end Integer scalar, end of the novel exon.
#' @param rstart Integer scalar, start of the downstream exon.
#' @param strand Factor scalar, strand of the novel exon (either "+" or "-").
#'
#' @return Character vector with transcript IDs.
#' @export
#' 
get_transcripts <- function(seqn, lend, start, end, rstart, strand, exons) {
  terminal <- 0
  ## if lend or rstart is NA, the exon is terminal
  if (!is.na(lend)) {
    ## find all exons that end in lend
    id1 <- mcols(subsetByOverlaps(exons,
                                  GRanges(seqnames = seqn,
                                          ranges = IRanges(lend, lend),
                                          strand = strand),
                                  type = "end",
                                  ignore.strand = FALSE))$transcript_id
  } else {
    terminal <- 1
  }
  if (!is.na(rstart)) {
    ## find all exons that start in rstart
    id2 <- mcols(subsetByOverlaps(exons,
                                  GRanges(seqnames = seqn,
                                          ranges = IRanges(rstart, rstart),
                                          strand = strand),
                                  type = "start",
                                  ignore.strand = FALSE))$transcript_id
  } else {
    terminal <- 2
  }
  if(terminal == 0) {   ## cassette exon
    ## make sure that the transcripts have no overlapping exon in the region
    ## where the new one will be located
    tr <- intersect(id1, id2)
    if (length(tr) == 0 ) NA
    
    olap_tr <- mcols(subsetByOverlaps(exons[mcols(exons)$transcript_id %in% tr],
                                      GRanges(seqnames = seqn,
                                              ranges = IRanges(lend + 1,
                                                               rstart - 1),
                                              strand = strand),
                                      ignore.strand = FALSE))$transcript_id
    ## all transcripts that do not have an exon in the region where the new exon
    ## will be located
    return(tr[!tr %in% olap_tr])
  } else if(terminal == 1) {   ## 5' terminal exon
    ## TODO: make sure that we only keep transcripts that start with the
    ## overlapping exons!!
    tr <- unique(id2)
    if (length(tr) == 0) NA
    
    olap_tr <- mcols(subsetByOverlaps(exons[mcols(exons)$transcript_id %in% tr],
                                      GRanges(seqnames = seqn,
                                              ranges = IRanges(start,
                                                               rstart - 1),
                                              strand = strand),
                                      ignore.strand = FALSE))$transcript_id
    return(tr[!tr %in% olap_tr])
  } else{  ## 3' terminal exon
    tr <- unique(id1)
    if (length(tr) == 0) NA
    ## remove all transcripts that have an exon in the connected intron or the
    ## novel exon
    olap_tr <- mcols(subsetByOverlaps(exons[mcols(exons)$transcript_id %in% tr],
                                      GRanges(seqnames = seqn,
                                              ranges = IRanges(lend + 1, end),
                                              strand = strand),
                                      ignore.strand = FALSE))$transcript_id
    return(tr[!tr %in% olap_tr])
  }
}


#' Create a new GTF entry for a novel exon
#'
#' Create GTF annotations for a novel exon, given a list of transcript IDs and
#' the GTF annotations of all exons. A random exon from each transcript is
#' copied and the start and end coordinates are exchanged with those of the
#' novel exon. The exon number and exon version are set to NA and the new exon
#' ID is the original exon id with a suffix that identifies the novel exon:
#' ID_me/exon_start:end
#'
#' @param tr_ids Character vector with transcrip ids as returned by
#'   [get_transcripts()]
#' @param exons GRanges object with exon annotations from a GTF file.
#' @param start Integer scalar, start of the novel exon.
#' @param end Integer scalar, end of the novel exon.
#'
#' @return GRanges object with the exon annotation of the novel exon and all
#'   exons from the transcripts in tr_ids.
#' @export
#' 
new_transcript <- function(start, end, tr_ids, exons) {
  type <- ifelse ((end - start + 1) < 28, "me", "exon")
  new_entries <- GRanges()
  
  tr_copy <- exons[which(mcols(exons)$transcript_id  %in% tr_ids)]
  mcols(tr_copy)$"transcript_id" <- paste0(mcols(tr_copy)$"transcript_id",
                                           "_", type, "_", start, ":", end)
  mcols(tr_copy)$"transcript_name" <- paste0(mcols(tr_copy)$"transcript_name",
                                             "_", type, "_", start, ":", end)
  
  copy <- tr_copy[match(paste0(tr_ids, "_", type, "_", start, ":", end),
                        mcols(tr_copy)$transcript_id )]
  ranges(copy) <- IRanges(start, end)  ## replace the ranges
  ## remove the exon id, number and version
  mcols(copy)[c("exon_number", "exon_version")] <- NA
  mcols(copy)$"exon_id" <- paste0(mcols(copy)$"exon_id", "_", type, "_",
                                  start, ":", end)
  
  return(c(tr_copy, copy))
}



#' Extend GTF file with predicted exons
#'
#' Add predicted novel exons to the correct transcripts in a GTF annotation.
#'
#' For each novel exon in `pred` (as returned by [find_novel_exons()], the list
#' of transcripts that contain its up- and downstream exons is determined. A
#' random exon from each transcript is copied and the start and end coordinates
#' are exchanged with those of the novel exon. The exon number and exon version
#' are set to NA and the new exon ID is the original exon id with a suffix that
#' identifies the novel exon: ID_me/exon_start:end
#'
#' @param gtf GTF annotations, either the path to the GTF file or a GRanges
#'   object
#' @param pred data.frame with predicted novel exons as returned by
#'   [find_novel_exons()]
#'
#' @return GRanges object with annotations from the GTF file and extended with
#'   the novel exons.
#' @importFrom rtracklayer import
#' @export
#' 
#' @examples 
#' gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery", 
#'                    mustWork = TRUE)
#'                    
#' ## Here we artificially create a data.frame with predicted exons. 
#' ## In general, the novel exons will be predicted with the function 
#' ## find_novel_exons()               
#' novel_exons <- data.frame(seqnames = c("19", "22"), 
#'                           start = c(47228064,41737092),
#'                           end = c(47228185, 41737150), strand = c("-", "+"),
#'                           lend = c(47226541, 41736141), 
#'                           rstart = c(47228589, 41738533))
#'                           
#' ## add the predicted exons to the GTF file                          
#' new_gtf <- extend_gtf(gtf, novel_exons)
#' 
#' ## save the new GTF to file with the export() function from rtracklayer
#' library(rtracklayer)
#' export(object = new_gtf, con = "extended_annotation.gtf")
extend_gtf <- function(gtf, pred) {
  if (is.character(gtf)) {
    gtf <- import(gtf)
  }
  exons <- gtf[mcols(gtf)$type == "exon", ]
  
  ## transcript IDs for each novel exon
  if (nrow(pred) == 1) {
    tr_ids <- list(get_transcripts(pred$seqnames, pred$lend, pred$start,
                                   pred$end, pred$rstart, pred$strand,exons))
  } else {
    ## TODO: speed up
    tr_ids <- mapply(get_transcripts, pred$seqnames, pred$lend, pred$start,
                     pred$end, pred$rstart, pred$strand,
                     MoreArgs = list(exons = exons))
  }
  ## remove the novel exons that do not have a transcript
  ## TODO: what to do with the removed predictions?
  pred <- pred[lengths(tr_ids) > 0, ]
  tr_ids <- tr_ids[lengths(tr_ids) > 0]
  
  novel_entries <- mapply(new_transcript, pred$start, pred$end, tr_ids,
                          MoreArgs = list(exons = exons))
  
  ## extend the original gtf with the novel entries
  c(gtf, do.call("c", novel_entries))
}
