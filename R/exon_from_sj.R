
#' Count mapped parts in a read from max-1 to 0 in the most downstream part
#' (from 5' to 3')
#'
#' @param x GRangesList object with mapped parts of reads
#'
#' @return data.frame with the range of each mapped part and the nr of each part
#'   from max-1 to 0 (most downstream part in the read)
#' @export
#'
downstream_count <- function(x) {
  tmp <- data.frame(x)
  tmp$nr <- nrow(tmp) - (1:nrow(tmp))
  tmp
}

#' Determine the coordinates of a novel exon at the end of a splice junction
#'
#' Given a set of reads that share a novel splice junction, determine the the
#' coordinates of a novel exon at the end of the splice junction. Only reads
#' with the splice junctions (the novel SJ and a downstream SJ) are considered.
#' If there are no reads with a second downstream SJ, take the range of the
#' longest mapped part as an approximation for the coordinates of the novel
#' exon (it probably is a terminal exon).
#'
#' @param r GAlignment object with reads that contain the novel splice
#'   junction
#' @param j_start integer start of the splice junction
#' @param j_end integer end of the splice junction
#' @param j_seqnames factor seqname of the splice junction
#' @param j_strand factor strand of the splice junction
#'
#' @return data.frame with the coordinates of the predicted exon(s)
#' @export
#'
identify_exon_end <- function(r, j_start, j_end, j_seqnames, j_strand) {
  r_ranges <- dplyr::bind_rows(lapply(r, downstream_count))

  sj_hit <- which(r_ranges$start == j_end + 1)
  sj_hit_downstream <- sj_hit[sj_hit %in% which(r_ranges$nr > 0)] + 1
  if (length(sj_hit_downstream) > 0) {
    ## If there are reads with a second junction, we take the coordinates from
    ## the mapped part.
    unique(data.frame(seqnames = j_seqnames,
                      lend = j_start - 1,
                      start = j_end + 1,
                      end = r_ranges$end[sj_hit_downstream - 1],
                      rstart = r_ranges$start[sj_hit_downstream],
                      strand = j_strand))
  } else {
    ## If there are no reads with a second junction, we take the longest mapped
    ## range.
    data.frame(seqnames = j_seqnames, lend = j_start - 1,
               start = j_end + 1,
               end = max(r_ranges$end[sj_hit]),
               rstart = NA,
               strand = j_strand)
  }
}



#' Count mapped parts in a read from 1 to max (from 5' to 3')
#'
#' @param x GRangesList object with mapped parts of reads
#'
#' @return data.frame with the range of each mapped part and the nr of each part
#'   from  1 to max (from 5' to 3')
#' @export
upstream_count <- function(x) {
  tmp <- data.frame(x)
  tmp$nr <- 1:nrow(tmp)
  tmp
}


#' Determine the coordinates of a novel exon at the start of a splice junction
#'
#' Given a set of reads that share a novel splice junction, determine the the
#' coordinates of a novel exon at the start of the splice junction. Only reads
#' with the splice junctions (the novel SJ and an upstream SJ) are considered.
#' If there are no reads with a second upstream SJ, take the range of the
#' longest mapped part as an approximation for the coordinates of the novel
#' exon (it probably is a terminal exon).
#' @param r GAlignment object with reads that contain the novel splice
#'   junction
#' @param j_start integer start of the splice junction
#' @param j_end integer end of the splice junction
#' @param j_seqnames factor seqname of the splice junction
#' @param j_strand factor strand of the splice junction
#'
#' @return data.frame with the coordinates of the predicted exon(s)
#' @export
#'
identify_exon_start <- function(r, j_start, j_end, j_seqnames, j_strand) {
  r_ranges <- dplyr::bind_rows(lapply(r, upstream_count))

  sj_hit <- which(r_ranges$end == j_start - 1)
  sj_hit_upstream <- sj_hit[r_ranges$nr[sj_hit] > 1]
  if (length(sj_hit_upstream) > 0) {
    ## if there are reads with an junction upstream of the SJ of interest
    unique(data.frame(seqnames = j_seqnames,
                      lend = r_ranges$end[sj_hit_upstream - 1],
                      start = r_ranges$start[sj_hit_upstream],
                      end = j_start - 1,
                      rstart = j_end + 1,
                      strand = j_strand))
  } else {
  ## If there are no reads with a second junction, we take the longest mapped
  ## range.
  data.frame(seqnames = j_seqnames,
             lend = NA,
             start = min(r_ranges$start[sj_hit]),
             end = j_start - 1,
             rstart = j_end + 1,
             strand = j_strand)
  }
}


#' Filter novel exon predictions of potential terminal exon.
#'
#' The exon predictions based on reads at the start or end of the novel splice
#' junction, are filtered and in case there are no reads with two splice
#' junctions (one of the coordinates is NA), determine if the novel exon could
#' be terminal. If yes, take the exon predictions from the corresponding end of
#' the splice junction.
#'
#' @param start_coords data.frame with exon predictions at the start of the
#'   novel SJ
#' @param end_coords data.frame with exon predictions at the end of the novel SJ
#' @param j data.frame with one row: the novel splice junction
#' @param txdb TxDb object, e.g. the txdb slot from the prepare_annotation()
#'   return object.
#' @param gtxdb GRanges: All genes from the TxDB object.
#' @param ebyTr GRanges: All exons from the TxDB object per transcript.
#'
#' @return data.frame with exon predictions, NULL if ambiguous or there are no
#'   supporting reads
#' @export
#'
filter_terminal_sj <- function(start_coords, end_coords, j,
                               txdb = txdb, gtxdb = gtxdb, ebyTr = ebyTr) {
  s_na <- is.na(start_coords$lend)
  e_na <- is.na(end_coords$rstart)
  if (any(!s_na)) {
    if (any(!e_na)) {
      dplyr::bind_rows(start_coords[!is.na(s_na), ],
                        end_coords[!is.na(e_na), ])
    } else {
      start_coords[!is.na(s_na), ]
    }
  } else if (any(!e_na)) {
      end_coords[!is.na(e_na), ]
  } else { ## terminal exon?
    terminal <- which_exon_terminal(GRanges(j), txdb = txdb, gtxdb = gtxdb,
                                    ebyTr = ebyTr)
    if (is.na(terminal)) {
      NULL
    } else if (terminal == "start") {
      end_coords
    } else if (terminal == "end") {
      start_coords
    }
  }
}

#' Identify the second splice junction of a novel exon
#'
#' This function takes a set of novel splice junctions (SJ) as input and returns
#' a data.frame with the location of the novel exon and the coordinates of the
#' neighbouring exons. There are three different cases: 1) Novel SJs that touch
#' an annotated exon with their start (5' end). 2) Novel SJs that touch
#' annotated exons with their end (3' end). 3) Novel SJs that touch and annoated
#' exon on both ends (5' and 3' end).
#'
#' The three cases can be illustrated as follows:
#' 1) Start touches annotation
#' AAAA            annotated exon
#'    J---J        novel junction
#'  xxx---xx----x  read
#'         X----X  function return
#'
#' 2) End touches annotation
#'            AAA annotated exon
#'       J----J   novel junction
#' xx---xx----xx  read
#'  X---X         function return
#'
#' 3) Both ends touch annotation
#'        AAAAA  annotated exon
#' AAAA          annotated exon
#'    J---J      novel junction
#'   nn---xxxxx  possible transcript with novel exon nn at the 5' of the SJ
#' xxxx---nn     possible transcript with novel exon nn at the 3' of the SJ
#'
#' In case 3, we search for reads with two novel splice junctions that support a
#' novel exon on either end of the SJ. If we do not find such reads, we check if
#' the novel exon could be terminal, i.e. the first or last exon in the
#' transcript. If yes, the novel exon does not have any touching annotated exon
#' on that end and thus we cannot determine the end coordinate of the novel
#' exon. As an approximation, we take the boundaries of the read with the
#' longest mapping to the novel exon.
#' @param junctions data.frame with novel splice junctions
#' @param reads GAlignment object with reads that contain novel splice
#'   junctions.
#' @param touching string defining which end of the novel splice junction
#'   touches an annotated exon.
#' @param txdb TxDb object, e.g. the txdb slot from the prepare_annotation()
#'   return object.
#' @param gtxdb GRanges: All genes from the TxDB object.
#' @param ebyTr GRanges: All exons from the TxDB object per transcript.
#'
#' @return data.frame with the coordinates of the novel exon. It has 6 columns:
#'   seqnames, lend, start, end, rstart and strand
#' @export
#'
get_second_sj <- function(junctions, reads, touching, txdb, gtxdb, ebyTr) {
  stopifnot(touching %in% c("start", "end", "both"))

  rs <- lapply(junctions$id, function(x) reads[mcols(reads)$which_label == x])
  r_mapped <- lapply(rs, function(x) grglist(x, order.as.in.query = FALSE,
                                             use.mcols = TRUE))
  ## remove all junctions without any supporting reads
  r_ind <- lengths(r_mapped) > 0
  if (!all(r_ind)) {
    junctions <- junctions[r_ind,]
    r_mapped <- r_mapped[r_ind]
  }
  if (touching == "both") {
    ## Case 3: The novel SJ touches annotated exons on both ends. We try to
    ## identify a novel exon on both the start and the end of the novel SJ
    full_coord_end <- mapply(identify_exon_end, r_mapped, junctions$start,
                             junctions$end, junctions$seqnames,
                             junctions$strand, SIMPLIFY = FALSE)
    full_coord_start <- mapply(identify_exon_start, r_mapped, junctions$start,
                               junctions$end, junctions$seqnames,
                               junctions$strand, SIMPLIFY = FALSE)
    ## apply functions are not working on GRanges anymore
    js_gr <- as.list(split(GRanges(junctions), 1:nrow(junctions)))

    ## We determine if any of the reads spanned two splice junctions. If not, we
    ## check if the novel exon could be terminal.
    full_coord <- mapply(filter_terminal_sj,
                         full_coord_start, full_coord_end, j = js_gr,
                         MoreArgs = list(txdb = txdb, gtxdb = gtxdb,
                                         ebyTr=ebyTr))
    return(dplyr::bind_rows(full_coord))
  } else if (touching == "start"){
    ## Case 1: The SJ touches an annotated exon on the start, so we check for
    ## reads that cover the novel SJ and the downstream exon.
    f <-  identify_exon_end
  } else {  ## "end"
    ## Case 2: The SJ touches an annotated exon on the end, so we check for
    ## reads that cover the novel SJ and the upstream exon.
    f <- identify_exon_start
  }
  full_coord <- mapply(f, r_mapped, junctions$start,
                       junctions$end, junctions$seqnames,
                       junctions$strand, SIMPLIFY = FALSE)
  dplyr::bind_rows(full_coord)
}



#' Identify a novel exon from a single novel splice junction
#'
#' Determine the size of the novel exon (unpaired internal splice junctions (SJ)
#' and terminal junctions) based on the touching anntated exons and the RNA-seq
#' reads. If a novel SJ touches either the start or end of an annotated exon and
#' we have supporting reads with two SJs (the novel SJ and an annotated one)
#' then we know the size of the novel exon and the coordinates of the
#' neighbouring exons. If a novel SJ touches annotated exons on both ends, then
#' the novel exon is either terminal or it overlaps with existing exons. In case
#' of a terminal exon, the function returns NA for the location of the
#' neighbouring exon.
#' @param reads GAlignment object with reads that contain novel splice
#'   junctions.
#' @param txdb TxDb object, e.g. the txdb slot from the prepare_annotation()
#'   return object.
#' @param gtxdb GRanges: All genes from the TxDB object.
#' @param ebyTr GRanges: All exons from the TxDB object per transcript.
#' @param df data.frame with novel splice junctions
#'
#' @return data.frame with the coordinates of the identifies novel exon from all
#'   splice junctions in df. It has 6 columns: seqnames, lend, start, end,
#'   rstart and strand
#'
#' @importFrom dplyr bind_rows
#'
#' @export
identify_exon_from_sj <- function(df, reads, txdb, gtxdb, ebyTr) {
  df$id <-  paste0(df$seqnames, ":", df$start, "-", df$end)

  ## split the junctions by type: "start", "end" or "both"
  js <- split(df, factor(df$touching))
  res <-  lapply(names(js), function(x) get_second_sj(js[[x]], reads, x, txdb,
                                                      gtxdb, ebyTr))
  dplyr::bind_rows(res)
}
