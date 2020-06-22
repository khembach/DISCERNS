#' DISCERNS
#'
#' The DISCERNS package predicts novel exons from RNA-seq data based on
#' the SJ.out.tab and the BAM file from STAR. The mapped splice-junctions are
#' compared with the genome annotation to discover novel exons. Existing genome
#' annotations in GTF format can be extended with the predicted exons.
#'
#' @docType package
#' @name DISCERNS
#' @author Katharina M. Hembach
#' @import Rcpp
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib DISCERNS, .registration = TRUE
NULL

