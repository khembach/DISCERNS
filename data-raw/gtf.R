library(rtracklayer)
library(dplyr)

base_dir <- "/home/Shared_sherborne/kathi/microexon_pipeline/simulation"
star_dir <- file.path(base_dir, "mapping/STAR/me_exon/outSJfilterOverhangMin6")
GTF <- file.path(base_dir, "reduced_GTF",
                 "GRCh37.85_chr19_22_reduced_me_exon.gtf")
SJFILE <- file.path(star_dir, "pass2_SJ.out.tab")
BAM <- file.path(star_dir, "pass2_Aligned.out_s.bam")

## list of removed exons
r_file <- file.path(base_dir, "reduced_GTF",
                    "removed_me_exon_unique_classified.txt")
dat <- read.table(r_file, header = TRUE, stringsAsFactors = FALSE)

## Select 6 exons: exons and microexons from the easy, medium and hard classes
## with more than one supporting read.
## We pick the first entry of each class
sel <- dat %>% group_by(type, class) %>%
  filter(count_reads > 1) %>%
  slice(1) %>%
  ungroup()

gtf <- import(GTF)
gtf <- gtf[gtf$gene_id %in% sel$gene_id]
usethis::use_data(gtf)
export(gtf, con = "inst/extdata/selected.gtf")


## Keep all SJs that are located within one of the selected genes
sj <- read.table(SJFILE)
## all SJs that overlap with one of the 6 genes
olap <- findOverlaps(GRanges(seqnames = sj$V1, ranges = IRanges(sj$V2, sj$V3),
                         strand = c("*", "+", "-")[sj$V4 + 1]),
                     gtf[gtf$type == "gene"],
                     ignore.strand = TRUE)
sj <- sj[unique(queryHits(olap)), ]
# usethis::use_data(sj)
write.table(sj, file = "inst/extdata/selected.SJ.out.tab", col.names = FALSE,
            row.names = FALSE, quote = FALSE)


## Select all reads from the BAM file that overlap any of the 6 selected exons
r <- paste0("'", paste0(as.character(GRanges(sel), ignore.strand = TRUE),
                        collapse = "' '"), "'")
## the read ids of all overlapping reads
cmd <- paste0("samtools view ", BAM, " ", r, " | cut -f1")
ids <- fread(cmd = cmd)
## filter all paired-end reads from the BAM file
## from https://www.biostars.org/p/105714/
cmd <- paste0("samtools view -h ", BAM, " | egrep '^@|^#|",
               paste0(unique(unlist(ids)), collapse = "|"),
               "' > tmp.sam")
system(cmd)
system("samtools view -b tmp.sam > inst/extdata/selected.bam")
system("rm tmp.sam")
system("samtools index inst/extdata/selected.bam")
