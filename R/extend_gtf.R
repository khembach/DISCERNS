#' Create a new GTF entry for a novel terminal exon
#'
#' This function creates new GTF entries for a novel terminal exon and all
#' transcripts in which it is located. As input, the function takes a list of
#' transcript IDs and the GTF annotations of all exons from the transcripts. A
#' random exon from each transcript is copied and the start and end coordinates
#' are exchanged with those of the novel exon. The ranges of each transcript are
#' adjusted to include the novel exon. The `exon_number` and `exon_version` are
#' set to NA. The new `exon_id`, `transcript_id` and `transcript_name` are the
#' original values with a suffix that identifies the novel exon by its ID: e.g.
#' `exon_id`_ID or `transcript_name`_ID.
#'
#' @param start Integer scalar. Start of the novel exon.
#' @param end Integer scalar. End of the novel exon.
#' @param pred_id Integer scalar. ID of the predicted exon.
#' @param tr_ids Character vector. Transcript IDs of all transcripts that
#'   contain the novel exon.
#' @param exon_copy GRanges object. All exons from the transcripts that contain
#'   the novel exon.
#' @param tran GRanges object with all transcript annotations from the organism.
#' 
#' @return GRanges object with the exon and transcript annotation entries of the
#'   novel exon and all exons from the transcripts in tr_ids.
#' @export
#'
create_tr_entry <- function(start, end, pred_id, tr_ids, exon_copy,
                            tran){
  ## create entry for new exon
  new_exon <- exon_copy[match(tr_ids, mcols(exon_copy)$transcript_id)]
  ranges(new_exon) <- IRanges(start, end)
  ## remove the exon number and version and update the exon_id
  mcols(new_exon)[c("exon_number", "exon_version")] <- NA
  mcols(new_exon)$"exon_id" <- paste0(mcols(new_exon)$"exon_id", "_", pred_id)

  tran <- tran[which(mcols(tran)$transcript_id %in% tr_ids)]
  ## compute the range of each transcript after adding the novel exon
  new_tr <- c(exon_copy, new_exon)
  tr_ranges <- lapply(split(new_tr, mcols(new_tr)$transcript_id), range)
  tr_ranges <- ranges(unlist(GRangesList(tr_ranges)))
  ## update the transcript ranges
  ranges(tran[match(names(tr_ranges), mcols(tran)$transcript_id)]) <- tr_ranges

  ## change the transcript_id and transcript_name of the trs and all their exons
  new_tr <- c(new_tr, tran)
  mcols(new_tr)$"transcript_id" <- paste0(mcols(new_tr)$"transcript_id",
                                          "_", pred_id)
  mcols(new_tr)$"transcript_name" <- paste0(mcols(new_tr)$"transcript_name",
                                            "_", pred_id)
  new_tr
}



#' Create new GTF annotation for a novel left terminal exon
#' 
#' Given a novel exon, the function identifies the IDs of all transcripts that
#' contain the downstream neighbouring exon (exon start is `rstart`) of the
#' novel exon. Each of the transcripts is copied and an entry for the the novel
#' exon is created. The `exon_number` and `exon_version` are set to NA. The new
#' `exon_id`, `transcript_id` and `transcript_name` are the original values with
#' a suffix that identifies the novel exon by its ID: e.g. `exon_id`_ID or
#' `transcript_name`_ID.
#'
#' @param seqn Character string or factor. seqname of the novel exon.
#' @param start Integer scalar. Start of the novel exon.
#' @param end Integer scalar. End of the novel exon.
#' @param rstart Integer scalar. Start of the downstream exon.
#' @param strand Character string or factor. Strand of the novel exon (either
#'   "+" or "-").
#' @param pred_id Integer scalar. ID of the novel exon.
#' @param exons GRanges object. Exon annotations from a GTF file.
#' @param tran GRanges object. Transcript annotations from a GTF file.
#' 
#' @return GRanges object with GTF annotations of the novel exon and all its
#'   transcripts
#'   
#' @export
#' 
new_left_terminal_transcript <- function(seqn, start, end, rstart, strand, 
                                         pred_id, exons, tran) {
  ## find all exons that start in rstart
  tr <- unique(mcols(subsetByOverlaps(exons, 
                                      GRanges(seqnames = seqn,
                                              ranges = IRanges(rstart, rstart),
                                              strand = strand),
                                      type = "start",
                                      ignore.strand = FALSE))$transcript_id)
  if (length(tr) > 0) {
    ## Are there any transcrips with exons upstream of rstart?
    exons <- exons[mcols(exons)$transcript_id %in% tr]
    olap_tr <- unique(exons[start(exons) < rstart]$transcript_id)
    if (length(olap_tr) == length(tr)) {
      exon_copy <- exons[which(mcols(exons)$transcript_id %in% tr)]
      ## remove all exons upstream of the neighbouring exon
      exon_copy <- exon_copy[!start(exon_copy) < rstart]
      create_tr_entry(start, end, pred_id, tr, exon_copy, tran)
    } else {
      tr <- tr[!tr %in% olap_tr]
      exon_copy <- exons[which(mcols(exons)$transcript_id %in% tr)]
      ## all transcripts start with the neighbouring exon
      create_tr_entry(start, end, pred_id, tr, exon_copy, tran)
    }
  }
}



#' Create new GTF annotation for a novel right terminal exon
#' 
#' Given a novel exon, the function identifies the IDs of all transcripts that
#' contain the upstream neighbouring exon (exon end is `lend`) of the
#' novel exon. Each of the transcripts is copied and an entry for the the novel
#' exon is created. The `exon_number` and `exon_version` are set to NA. The new
#' `exon_id`, `transcript_id` and `transcript_name` are the original values with
#' a suffix that identifies the novel exon by its ID: e.g. `exon_id`_ID or
#' `transcript_name`_ID.
#'
#' @param seqn Character string or factor. seqname of the novel exon.
#' @param lend Integer scalar. End of the upstream exon.
#' @param start Integer scalar. Start of the novel exon.
#' @param end Integer scalar. End of the novel exon.
#' @param strand Character string or factor. Strand of the novel exon (either
#'   "+" or "-").
#' @param pred_id Integer scalar. ID of the novel exon.
#' @param exons GRanges object. Exon annotations from a GTF file.
#' @param tran GRanges object. Transcript annotations from a GTF file.
#' 
#' @return GRanges object with GTF annotations of the novel exon and all its
#'   transcripts
#'   
#' @export
#'
new_right_terminal_transcript <- function(seqn, lend, start, end, strand, 
                                          pred_id, exons, tran) {
  ## find all exons that end in lend
  tr <- unique(mcols(subsetByOverlaps(exons, 
                                      GRanges(seqnames = seqn,
                                              ranges = IRanges(lend, lend),
                                              strand = strand),
                                      type = "end",
                                      ignore.strand = FALSE))$transcript_id)
  if (length(tr) > 0) {
    ## Are there any transcrips with exons downtream of lend?
    exons <- exons[mcols(exons)$transcript_id %in% tr]
    olap_tr <- unique(exons[start(exons) > lend]$transcript_id)
    
    if (length(olap_tr) == length(tr)) {  ## all tr have overlapping exons
      exon_copy <- exons[which(mcols(exons)$transcript_id %in% tr)]
      ## remove all downstream exons of the neighbouring exon
      exon_copy <- exon_copy[!end(exon_copy) > lend]
      create_tr_entry(start, end, pred_id, tr, exon_copy, tran)
    } else {
      tr <- tr[!tr %in% olap_tr]
      exon_copy <- exons[which(mcols(exons)$transcript_id %in% tr)]
      ## all transcripts end in the neighbouring exon
      create_tr_entry(start, end, pred_id, tr, exon_copy, tran)
    }
  }
}


#' Create new GTF annotation for a novel internal exon
#' 
#' Given a novel exon, the function identifies the IDs of all transcripts that
#' contain the upstream and downstream neighbouring exons (exon end is `lend`
#' and exon start is `rstart`) of the novel exon. Each of the transcripts is
#' copied and an entry for the the novel exon is created. The `exon_number` and
#' `exon_version` are set to NA. The new `exon_id`, `transcript_id` and
#' `transcript_name` are the original values with a suffix that identifies the
#' novel exon by its ID: e.g. `exon_id`_ID or
#' `transcript_name`_ID.
#'
#' @param seqn Character string or factor. seqname of the novel exon.
#' @param lend Integer scalar. End of the upstream exon.
#' @param start Integer scalar. Start of the novel exon.
#' @param end Integer scalar. End of the novel exon.
#' @param rstart Integer scalar. Start of the downstream exon.
#' @param strand Character string or factor. Strand of the novel exon (either
#'   "+" or "-").
#' @param pred_id Integer scalar. ID of the novel exon.
#' @param exons GRanges object. Exon annotations from a GTF file.
#' @param tran GRanges object. Transcript annotations from a GTF file.
#' 
#' @return GRanges object with GTF annotations of the novel exon and all its
#'   transcripts
#'   
#' @export
#' 
new_internal_transcript <- function(seqn, lend, start, end, rstart, strand, 
                                    pred_id, exons, tran) {
  ## ids of all transcripts with exons that end in lend and start in rstart
  id1 <- mcols(subsetByOverlaps(exons,
                                GRanges(seqnames = seqn,
                                        ranges = IRanges(lend, lend),
                                        strand = strand),
                                type = "end",
                                ignore.strand = FALSE))$transcript_id
  id2 <- mcols(subsetByOverlaps(exons,
                                GRanges(seqnames = seqn,
                                        ranges = IRanges(rstart, rstart),
                                        strand = strand),
                                type = "start",
                                ignore.strand = FALSE))$transcript_id
  tr <- intersect(id1, id2)
  if (length(tr) > 0 ) {
    ## we need a transcript that has an intron from lend to rstart
    ## Are there any transcripts with an exon in the intron?
    olap_tr <- unique(mcols(
      subsetByOverlaps(exons[mcols(exons)$transcript_id %in% tr],
                       GRanges(seqnames = seqn, 
                               ranges = IRanges(lend + 1, rstart - 1),
                               strand = strand),
                       ignore.strand = FALSE)
      )$transcript_id)
    
    if (length(olap_tr) == length(tr)) { 
      ## all transcripts have an exon in the intron
      exon_copy <- exons[which(mcols(exons)$transcript_id %in% tr)]
      ## all exons inside the intron
      olap_exons <- subsetByOverlaps(exon_copy,
                                     GRanges(seqnames = seqn,
                                             ranges = IRanges(lend + 1,
                                                              rstart - 1),
                                             strand = strand),
                                     ignore.strand = FALSE)
      ## take a random exon per transcript and replace its coordinates
      new_exon <- olap_exons[match(tr, mcols(olap_exons)$transcript_id)]
      ranges(new_exon) <- IRanges(start, end)
      ## remove the exon number and version and update the exon_id
      mcols(new_exon)[c("exon_number", "exon_version")] <- NA
      mcols(new_exon)$"exon_id" <- paste0(mcols(new_exon)$"exon_id", "_", 
                                          pred_id)
      ## create new transcript entries and change the transcript_id and
      ## transcript_name of the transcripts and all their exons
      new_tr <- c(exon_copy, new_exon, 
                  tran[which(mcols(tran)$transcript_id %in% tr)])
      mcols(new_tr)$"transcript_id" <- paste0(mcols(new_tr)$"transcript_id",
                                              "_", pred_id)
      mcols(new_tr)$"transcript_name" <- paste0(mcols(new_tr)$"transcript_name",
                                                "_", pred_id)
      new_tr 
    } else { ## all transcripts with no exon in the intron
      tr <- tr[!tr %in% olap_tr]
      exon_copy <- exons[which(mcols(exons)$transcript_id %in% tr)]
      ## take a random exon per transcript
      new_exon <- exon_copy[match(tr, mcols(exon_copy)$transcript_id)]
      ranges(new_exon) <- IRanges(start, end)
      ## remove the exon number and version and update the exon_id
      mcols(new_exon)[c("exon_number", "exon_version")] <- NA
      mcols(new_exon)$"exon_id" <- paste0(mcols(new_exon)$"exon_id", "_", 
                                          pred_id)
      ## create new transcript entries and change the transcript_id and
      ## transcript_name of the transcripts and all their exons
      new_tr <- c(exon_copy, new_exon, 
                  tran[which(mcols(tran)$transcript_id %in% tr)])
      mcols(new_tr)$"transcript_id" <- paste0(mcols(new_tr)$"transcript_id",
                                              "_", pred_id)
      mcols(new_tr)$"transcript_name" <- paste0(mcols(new_tr)$"transcript_name",
                                                "_", pred_id)
      new_tr 
    }
  }
}





#' Extend GTF file with predicted exons
#'
#' Add predicted novel exons to the correct transcripts in a GTF annotation.
#'
#' For each novel exon in `pred` (as returned by [find_novel_exons()]), the list
#' of transcripts that contain its up- and/or downstream exons is determined. A
#' random exon from each transcript is copied and the start and end coordinates
#' are exchanged with those of the novel exon. The `exon_number` and
#' `exon_version` are set to NA. The new `exon_id`, `transcript_id` and
#' `transcript_name` are the original values with a suffix that identifies the
#' novel exon by its ID: e.g. `exon_id`_ID or `transcript_name`_ID.
#'
#' @param gtf GTF annotations, either the path to the GTF file as a character
#'   string or a GRanges object.
#' @param pred data.frame with predicted novel exons as returned by
#'   [find_novel_exons()].
#' @param cores Integer scalar. Number of cores to use. Default 1.
#'
#' @return GRanges object with annotations from the GTF file and extended with
#'   the novel xons.
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
#'                           rstart = c(47228589, 41738533), ID = c(1, 2))
#'
#' ## add the predicted exons to the GTF file
#' new_gtf <- extend_gtf(gtf, novel_exons)
#'
#' ## save the new GTF to file with the export() function from rtracklayer
#' library(rtracklayer)
#' export(object = new_gtf, con = "extended_annotation.gtf")
extend_gtf <- function(gtf, pred, cores = 1) {
  if (is.character(gtf)) {
    gtf <- import(gtf)
  }
  exons <- gtf[mcols(gtf)$type == "exon", ]
  tran <- gtf[mcols(gtf)$type == "transcript", ]

  ## Seperate the cassette exons from the terminal exon, because we need to
  ## adjust the transcript range for terminal exons but not for internal exons.
  type <- ifelse(is.na(pred$lend), "left_terminal", 
                 ifelse(is.na(pred$rstart), "right_terminal", "internal"))
  pred_split <- split(pred, type)

  if ("left_terminal" %in% names(pred_split)) {
    p <- pred_split[["left_terminal"]]
    novel_entries <- mcmapply(new_left_terminal_transcript, p$seqnames, p$start, 
                              p$end, p$rstart, p$strand, p$ID,
                              MoreArgs = list(exons = exons, tran = tran), 
                              mc.cores = cores, SIMPLIFY = FALSE)
    gtf <- c(gtf, do.call("c", novel_entries[!sapply(novel_entries, is.null)]))
    
  }
  if ("right_terminal" %in% names(pred_split)) {
    p <- pred_split[["right_terminal"]]
    novel_entries <- mcmapply(new_right_terminal_transcript, p$seqnames, p$lend,
                              p$start, p$end, p$strand, p$ID,
                              MoreArgs = list(exons = exons, tran = tran), 
                              mc.cores = cores, SIMPLIFY = FALSE)
    gtf <- c(gtf, do.call("c", novel_entries[!sapply(novel_entries, is.null)]))
  }
  if ("internal" %in% names(pred_split)) {
    p <- pred_split[["internal"]]
    novel_entries <- mcmapply(new_internal_transcript, p$seqnames, p$lend,
                              p$start, p$end, p$rstart, p$strand, p$ID,
                              MoreArgs = list(exons = exons, tran = tran), 
                              mc.cores = cores, SIMPLIFY = FALSE)
    gtf <- c(gtf, do.call("c", novel_entries[!sapply(novel_entries, is.null)]))
  }
  gtf
}