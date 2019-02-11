# exondiscovery

The `exondiscovery` R package contains functions for the discovery of novel exons in RNA-seq data. The package requires read alignments from STAR because it is based on the `SJ.out.tab` file. Additionally, the user has to input gene annotations in GTF format and a BAM file with the aligned RNA-seq reads (optional). The output is a table with the genomic coordinates of predicted exons. Each exon prediction includes the end coordinates of the upstream exon and the coordinates of the downstream exon. The GTF annotation can be supplemented with the predicted exons.
