#' Novel exon prediction
#'
#' Predict novel exons based on the SJ.out.tab file from STAR and/or a BAM file.
#' The function has three different prediction modes:
#'   1. Cassette exons
#'   2. Novel exons from a single novel splice junction
#'   3. Novel exons from reads/read pairs with two novel splice junctions.
#'   
#' The cassette exon prediction is always preformed and the second and third 
#' mode can be turned on/off with the paramters `single_sj` and `read_based`.
#'
#' The three different prediction modes are explained in more detail below:
#' 
#' Novel cassette exons are predicted from pairs of novel splice junctions (SJ)
#' that are located within an annotated intron and share the start and end
#' coordinates of the intron:
#' \preformatted{
#' X---------X   annotated intron
#' x---x         novel SJ
#'       x---x   novel SJ
#' X---NNN---X   predicted cassette exon (N)
#' }
#' First, the novel SJs are filtered: Only SJs that are located within an annotated
#' intron and that share their start or end coordinates with the intron are
#' retained. All introns that share both their start and end with a novel SJ are
#' tested for cassette exons. If the two novel SJs within an intron do not
#' overlap and are on the same strand, a novel cassette exon is predicted.
#' 
#' The second mode (parameter `single_sj`) is based on single novel SJs as input.
#' The novel SJs can be divided in three different cases: 
#'    1. Novel SJs that touch an annotated exon with their start (5' end). 
#'    2. Novel SJs that touch annotated exons with their end (3' end). 
#'    3. Novel SJs that touch an annoated exon on both ends (5' and 3' end).
#'
#' The three cases can be illustrated as follows:
#'    1. Start touches annotation
#'    \preformatted{
#' AAAA            annotated exon
#'    J---J        novel SJ
#'  xxx---xx----x  read
#'        NN----X  novel exon and second SJ
#'    }
#'    2. End touches annotation
#'    \preformatted{
#'            AAA annotated exon
#'       J----J   novel SJ
#' xx---xx----xx  read
#'  X---NN        novel exon and second SJ
#'    }
#'    3. Both ends touch annotation
#'    \preformatted{
#'        AAAAA  annotated exon
#' AAAA          annotated exon
#'    J---J      novel junction
#'   NN---xxxxx  possible transcript with novel exon NN at the 5' of the SJ
#' xxxx---NN     possible transcript with novel exon NN at the 3' of the SJ
#'    }
#' 
#' For case 1 and 2, the function reads the BAM file and takes all reads that
#' support the novel SJ and have a second SJ. The coordinates of the second SJ,
#' and thus the coordinates of the novel exon are determined from the reads. In
#' case 3, the function searches for reads with two novel SJs that support a
#' possible novel exon on either end of the SJ. If no reads are found, the
#' function checks if the novel exon could be terminal, i.e. the first or last
#' exon in a transcript. If yes, the novel exon is not connected to an annotated
#' exon at its start/end and thus the start/end coordinate of the novel exon
#' cannot be determined clearly. As an approximation, the function takes the
#' boundaries of the read with the longest mapping to the novel exon.

#' The third prediction mode (parameter `read_based`) only uses reads and the
#' annotation as input. Novel exons are predicted from novel combinations of
#' already annotated SJs. The SJ pairs are defined by a single read with two SJs
#' or by read pairs where each read spans one junction. Each of the SJ pairs are
#' compared with the annotated SJs per transcript and already annotated SJ pairs
#' are removed. 

#'For paired-end reads, there is an additional requirement: The distance between
#'the end of the first read in a pair and the start of the second read has to
#'has to be ```< 2*(readlength-minOverhang) + minIntronSize```. Here,
#'`readlength` is the length of the reads, `minOverhang` is the minmal required
#'read overhang over a SJ of the alignment tool and `minIntronSize` is the
#'minimal required intron length of the alignment tool. For example, paired-end
#'reads with a lenght of 101 nts and a minimal overhang of 6 and a minimal
#'intron length of 21 allow a distance of at most 211 nucleotides between the
#'two SJs: 2*(101-6) + 21 = 211. If the distance between the two SJs exceeds the
#'limit, it cannot be guaranteed that the junctions are connected to the same
#'exon.
#'
#'Novel exons are then predicted from novel combinations of the identified SJ
#'pairs. A novel exon is only predicted if it is contained within the boundaries
#'of an annotated gene. This prevents the prediction of false positive exons
#'because of wrongly mapped reads.
#'
#'@param sj_filename Path to SJ.out.tab file.
#'@param annotation List with exon and intron annotation as GRanges. Created
#'  with [prepare_annotation()].
#'@param min_unique Minimal number of reads required to map over a splice
#'  junction.
#'@param verbose Print messages about the progress.
#'@param gzipped A logical scalar. Is the input file gzipped? Default FALSE.
#'@param bam Path to BAM file.
#'@param single_sj  A logical scalar. Should single novel splice junctions be
#'  used for prediction? Requires a BAM file.
#'@param read_based  A logical scalar. Should novel exons be predicted from
#'  reads or read-pairs with two novel splice junctions? Requires a BAM file.
#'@param yield_size Integer scalar. Read the BAM file in chunks of this size.
#'@param lib_type Character scalar. Type of the sequencing library: either "SE"
#'  (single-end) or "PE" (paired-end); default "PE".
#'@param stranded Character scalar. Strand type of the sequencing protocol:
#'  "unstranded" for unstranded protocols; "forward" or "reverse" for stranded
#'  protocols. "forward" refers to protocols where the first read comes from the
#'  forward (sense) strand; default "reverse". See the [Salmon
#'  documentation](https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype)
#'  for an explanation of the different fragment library types.
#'
#'@return Data.frame with the coordinates of the identified novel exons. Each
#'  row in the data.frame is a predicted novel exon. The columns are the
#'  chromosome (`seqnames`), the end of the upstream exon (`lend`) in the
#'  transcript, the start and end of the novel exon (`start` and `end`), the
#'  start of the downstream exon (`rstart`) in the transcript and the strand
#'  (`strand`). The last three columns are the number of reads supporting each
#'  of the two splice junctions that define the novel exon: `unique_left` is the
#'  number of reads supporting the SJ from `lend` to `start` and `unique_right`
#'  is the number of supporting reads for the SJ from `end` to `rstart`.
#'  `min_reads` is the minimum of the two.
#'
#'@importFrom data.table fread
#'@import GenomicRanges
#'@import IRanges
#'@importFrom S4Vectors queryHits subjectHits
#'@importFrom GenomicFeatures genes exonsBy
#'@importFrom dplyr full_join mutate
#'
#'@export
#'
#'@examples
#' sj <- system.file("extdata", "selected.SJ.out.tab",
#'                   package = "exondiscovery", mustWork = TRUE)
#' bam <- system.file("extdata", "selected.bam",
#'                    package = "exondiscovery", mustWork = TRUE)
#'
#' ## prepare annotation
#' gtf <- system.file("extdata", "selected.gtf", package = "exondiscovery",
#'                    mustWork = TRUE)
#' anno <- prepare_annotation(gtf)
#'
#' ## predict novel exons
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam)
#'
#' ## predict only cassette exons
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  single_sj = FALSE, read_based = FALSE)
#'
#' ## only predict exons from novel splice junctions
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, read_based = FALSE)
#'
#' ## Only consider novel splice junctions with at least ten supporting reads
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 10,
#'                  bam = bam)
#'
#' ## turn verbose off
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, verbose = FALSE)
#'
#' ## increase chunk size for BAM file reading if lots of memory is available
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, yield_size = 1000000)
#'
#' ## Stranded single-end reads from the forward strand (sense)
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, lib_type = "SE", stranded = "forward")
#'
#' ## Stranded paired-end reads where the first read comes from the forward strand
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, lib_type = "PE", stranded = "forward")
#'
#' ## Stranded paired-end reads where the first read comes from the reverse
#' ## strand, e.g., Illumina TruSeq stranded mRNA protocol. This is the default.
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, lib_type = "PE", stranded = "reverse")
#' 
#' ## Unstranded single-end reads
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, lib_type = "SE", stranded = "unstranded")
#' 
find_novel_exons <- function(sj_filename, annotation, min_unique = 1,
                             gzipped = FALSE, verbose = TRUE, bam,
                             single_sj = TRUE, read_based = TRUE, 
                             yield_size = 200000, 
                             lib_type = "PE", stranded = "reverse") {
  
  if (!lib_type %in% c("SE", "PE")) {
    stop('Parameter lib_type has to be either "SE" or "PE"')
  }
  if (!stranded %in% c("unstranded", "forward", "reverse")) {
    stop('Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  }
  
  exons <- annotation[["exons"]]
  introns <- annotation[["introns"]]
  
  if (verbose) message("Reading splice junctions from file")
  
  if (gzipped) {
    sj <- fread(paste0("zcat ", sj_filename))
  } else {
    sj <- fread(sj_filename)
  }
  colnames(sj) <- c("seqnames", "start", "end", "strand", "motif", "annotated",
                    "unique", "mutimapping", "maxoverhang")
  
  #strand: (0: undefined, 1: +, 2: -)
  sj$strand <- c("*", "+", "-")[sj$strand + 1]
  
  ### TODO: create functions for the different prediction "types"
  ## 1) casette exons 2) predictions from single novel SJ
  ## 3) predictions from reads or read pairs with 2 SJ
  ## let the user decide which ones to run
  
  ## ============= junction prefiltering ============= #
  if (verbose) message("Prefiltering splice junctions")
  
  sj_gr <- GRanges(sj)
  
  ## filter out all annotated junctions
  sj_ann <- subsetByOverlaps(sj_gr, introns, type = "equal")
  sj_unann <- sj_gr[!(sj_gr %in% sj_ann)]
  
  ## filter out all junctions with less than min_unique reads
  sj_unann <- sj_unann[mcols(sj_unann)$unique >= min_unique]
  
  ## filter out all novel junctions that do not touch an annotated exon on the
  ## same strand
  sj_unann <- sj_unann[(start(sj_unann) - 1) %in% end(exons) |
                         (end(sj_unann) + 1) %in% start(exons), ]
  
  ## ======== Cassette exon prediction ======== #
  
  if (verbose) message("Predicting cassette exons")
  ce <- predict_cassette_exon(sj_unann, introns)
  novel_exons <- ce[["ne"]]
  if (length(ce[["sj"]]) > 0) {
    sj_unann <- ce[["sj"]]
  } else {
    single_sj <- FALSE
    if (verbose) message("No novel SJs left, skipping single SJ step.")
  }
  
  ## ==== Find novel exon coordinates from single novel splice junctions ======
  if (single_sj) {
    if (verbose)
      message("Find exon coordinates of exons adjacent to novel splice junctions")
    if (missing(bam)) {
      stop("Please specify a BAM file with parameter bam.")
    }
    reads <- import_novel_sj_reads(bam, sj_unann)
    
    ebyTr <- exonsBy(annotation[["txdb"]], by = "tx", use.names = TRUE)
    gtxdb <- genes(annotation[["txdb"]])
    
    ## convert GRanges to data.frame otherwise we cannot use apply functions
    sj_unann <- data.frame(sj_unann, stringsAsFactors = TRUE)
    
    ## annotate which end of a splice junction is touching an annotated exon
    touching <- mapply(sj_touching_exon, sj_unann$seqnames, sj_unann$start,
                       sj_unann$end, sj_unann$strand,
                       MoreArgs = list(exons = exons))
    sj_unann$touching <- touching
    ## remove the ambiguous sj that touch different genes on both ends
    sj_unann <- sj_unann[!is.na(touching), ]
    
    res <- identify_exon_from_sj(sj_unann, reads = reads,
                                 txdb = annotation[["txdb"]],
                                 gtxdb = gtxdb, ebyTr = ebyTr)
    
    novel_exons <- rbind(novel_exons, res)
    
  }
  
  if (read_based) {
    ## ============ Predict novel exons from reads with 2 junctions ==============
    if (verbose) message("Predict novel exons from reads with 2 junctions")
    
    if (missing(bam)) {
      stop("Please specify a BAM file with parameter bam.")
    }
    junc_reads <- filter_junction_reads(bam, yield_size = yield_size, 
                                        lib_type = lib_type, 
                                        stranded = stranded)
    read_pred <- predict_jr_exon(junc_reads, annotation)
    
    if (exists("novel_exons", inherits = FALSE)) {
      ## make sure the factors have the same levels
      combined <- union(levels(read_pred$strand), levels(novel_exons$strand))
      ## predictions from both the reads and the SJ.out.tab
      novel_exons <- full_join(mutate(read_pred,
                                      strand = factor(strand, 
                                                      levels = combined)),
                               mutate(novel_exons,
                                      strand = factor(strand,
                                                      levels = combined)),
                               by = c("seqnames", "start", "end", "strand", 
                                      "lend", "rstart"))
    } else {
      novel_exons <- read_pred
    }
    
    ## ======= Predict novel exons from read pairs with each one junction =======
    if (verbose) message("Predict novel exons from read pairs with each 1 junction")
    
    read_pair_pred <- predict_jrp_exon(junc_reads, annotation)
    combined <- union(levels(read_pair_pred$strand), levels(novel_exons$strand))
    novel_exons <- full_join(mutate(read_pair_pred,
                                    strand = factor(strand, levels = combined)),
                             mutate(novel_exons,
                                    strand = factor(strand, levels = combined)),
                             by = c("seqnames", "start", "end", "strand", 
                                    "lend", "rstart"))
  }
  
  ## ============ Compute minimal junction read coverage ===========
  if(verbose) message("Computing minimal junction read coverage")
  
  ##  Add columns with the number of reads supporting the left and right splice
  ##  junction and the minimum of both
  key_sj <- with(sj, paste(seqnames, start - 1, end + 1, strand, sep = "."))
  key_ne_l <- paste(novel_exons$seqnames, novel_exons$lend, novel_exons$start,
                    novel_exons$strand, sep = ".")
  key_ne_r <- paste(novel_exons$seqnames, novel_exons$end, novel_exons$rstart,
                    novel_exons$strand, sep = ".")
  
  novel_exons$unique_left <- sj$unique[match(key_ne_l, key_sj)]
  novel_exons$unique_right <- sj$unique[match(key_ne_r, key_sj)]
  
  ## if the exon is a cassette exons, but one of the junctions is not in
  #SJ.out.tab (count = NA), # remove the prediction, because it is most likely
  #wrong (this only happens if we predict exons from reads) novel_exons <-
  #novel_exons %>% filter(!(!is.na(lend) & is.na(unique_left))) %>%
  #filter(!(!is.na(rstart) & is.na(unique_right)))
  
  ## take the minimum read coverage of both junctions
  novel_exons$min_reads <- pmin(novel_exons$unique_left,
                                novel_exons$unique_right, na.rm = TRUE)
  novel_exons$ID <- 1:nrow(novel_exons)
  novel_exons
}


