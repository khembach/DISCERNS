% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exon_from_sj.R
\name{filter_terminal_sj}
\alias{filter_terminal_sj}
\title{Filter novel exon predictions of potential terminal exon}
\usage{
filter_terminal_sj(
  start_coords,
  end_coords,
  j,
  txdb = txdb,
  gtxdb = gtxdb,
  ebyTr = ebyTr
)
}
\arguments{
\item{start_coords}{data.frame with exon predictions at the start of the
novel SJ.}

\item{end_coords}{data.frame with exon predictions at the end of the novel SJ.}

\item{j}{data.frame with one row: the novel splice junction.}

\item{txdb}{TxDb object, e.g. the "txdb" slot from the \code{\link[=prepare_annotation]{prepare_annotation()}}
return object.}

\item{gtxdb}{GRanges object. All genes from the txdb parameter, e.g. obtained
with \code{GenomicFeatures::genes(txdb)}.}

\item{ebyTr}{GRangesList object. All exons per transcript of the txdb
parameter, e.g. obtained with \code{GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)}.}
}
\value{
data.frame with exon predictions, NULL if ambiguous or there are no
supporting reads
}
\description{
The exon predictions based on reads at the start or end of the novel splice
junction are filtered and in case there are no reads with two splice
junctions (one of the coordinates is NA), determine if the novel exon could
be terminal. If yes, take the exon predictions from the corresponding end of
the splice junction.
}
