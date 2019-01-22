#' Determine which side of the splice junction touches an exon
#'
#' Find all exons with the same strand as the splice junction and that touch it.
#' If a sj touches two exons, the exons must come from the same gene.
#'
#' @param sj A splice junction as a single GRange
#' @param exons Exon annotation as GRanges
#'
#' @return A string that defines which side of the splice junction is touching
#'   an exon: "both", "start", "end", NA if the splic jucntion does not touch an
#'   exon or if the exons come from different genes.
#' @export
#'
sj_touching_exon <- function(sj, exons){
  s <- which((start(sj)-1) == end(exons) )
  s <- exons[s,]   ## all touching exons at start of sj
  s <- s[strand(s) == strand(sj),] ## all touching exons with same strand

  e <- which((end(sj)+1) == start(exons))
  e <- exons[e]   ## all touching exons at end of sj
  e <- e[strand(e) == strand(sj),]

  if( length(s) >0 ){
    if(length(e) >0){
      ## pair of touching exons from the same gene
      if(any( mcols(s)$gene_id %in% mcols(e)$gene_id )) {
        "both"
      } else NA
    } else "start"
  } else if(length(e) >0){
    "end"
  } else NA
}


#' Match novel splice junctions within an intron
#'
#' @param s Junctions from the start of the intron as GRanges
#' @param e Junctions from the end of the intron as GRanges
#'
#' @return A vector with the seqnames, end of the preceding exon, start of the
#'   novel exon, end of the novel exon, start of the consecutive exon, and the
#'   strand.
#' @export
#'
match_sj_in_intron <- function(s, e){
  if (all( length(s) ==1, length(e)==1,
           isDisjoint(c(s, e), ignore.strand=FALSE))){
    return(c(as.vector(seqnames(s)), start(s)-1, end(s)+1,
               start(e)-1, end(e)+1, as.vector(strand(s)) ) )
  } ## TODO: more than two novel splice junctions per intron
}
