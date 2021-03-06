% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_novel_sj_reads.R
\name{import_novel_sj_reads}
\alias{import_novel_sj_reads}
\title{Import reads with novel splice junctions}
\usage{
import_novel_sj_reads(bam, sj_unann)
}
\arguments{
\item{bam}{Character string. The path to the BAM file.}

\item{sj_unann}{GRanges object with unannotated splice junctions.}
}
\value{
GAlignments object with all reads from the BAM file that contain a
novel splice junctions.
}
\description{
This functions takes a BAM file and a GRanges object with unannotated splice
junctions as input and returns a GAlignments object with all reads from the
BAM file that contain any of the novel splice junctions.
}
