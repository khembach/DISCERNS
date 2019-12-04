#' Novel exon prediction
#'
#' Predict novel exons based on the SJ.out.tab file from STAR and/or a BAM file.
#' The function has three different prediction modes:
#'   1. Cassette exons
#'   2. Novel exons from a single novel splice junction
#'   3. Novel exons from reads/read pairs with two novel splice junctions.
#'   
#' The cassette exon prediction is always performed and the second and third 
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
#'the end of the first read in a pair and the start of the second read has to be
#'```< 2*(readlength-minOverhang) + minIntronSize```. Here, `readlength` is the
#'length of the reads, `minOverhang` is the minmal required read overhang over a
#'SJ of the alignment tool and `minIntronSize` is the minimal required intron
#'length of the alignment tool. For example, paired-end reads with a lenght of
#'101 nts and a minimal overhang of 6 and a minimal intron length of 21 allow a
#'distance of at most 211 nucleotides between the two SJs: 2*(101-6) + 21 = 211.
#'If the distance between the two SJs exceeds the limit, it cannot be guaranteed
#'that the junctions are connected to the same exon.
#'
#'Novel exons are then predicted from novel combinations of the identified SJ
#'pairs. A novel exon is only predicted if it is contained within the boundaries
#'of an annotated gene. This prevents the prediction of false positive exons
#'because of wrongly mapped reads.
#'
#'@param sj_filename Path to SJ.out.tab file.
#'@param annotation List with exon and intron annotation as GRanges. Created
#'  with [prepare_annotation()].
#'@param min_unique Integer scalar. Minimal number of reads required to map over
#'  a splice junction. Novel splice junctions with less reads are excluded from
#'  the analysis. All predicted novel exons with less than `min_unique` reads
#'  are discared.
#'@param verbose Print messages about the progress.
#'@param gzipped A logical scalar. Is the input file gzipped? Default FALSE.
#'@param bam Character string. The path to the BAM file.
#'@param single_sj  A logical scalar. Should single novel splice junctions be
#'  used for prediction? Requires a BAM file.
#'@param read_based  A logical scalar. Should novel exons be predicted from
#'  reads or read-pairs with two novel splice junctions? Requires a BAM file.
#'@param yield_size Integer scalar. Read the BAM file in chunks of this size.
#'@param lib_type Character string. Type of the sequencing library: either "SE"
#'  (single-end) or "PE" (paired-end). Default "PE".
#'@param stranded Character string. Strand type of the sequencing protocol:
#'  "unstranded" for unstranded protocols; "forward" or "reverse" for stranded
#'  protocols. In a "forward" protocol, the first read in a pair (or sinlge-end
#'  reads) comes from the forward (sense) strand and in a "reverse" protocol
#'  from the reverse (antisense) strand. See the [Salmon
#'  documentation](https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype)
#'   for an explanation of the different fragment library types. Default
#'  "reverse".
#' @param read_length Integer scalar. Length of your reads in bps. Default 101.
#' @param overhang_min Integer scalar. Minimum overhang length for splice
#'   junctions on both sides as defined by the `--outSJfilterOverhangMin`
#'   parameter of STAR. Use the minimum of the values for canonical splice junctions
#'   (value (2) to (4)). You do not have to set this parameter if you used the
#'   default values from STAR. Default 12.
#' @param min_intron_size Integer scalar. Minimum intron size
#'   (`--alignIntronMin` parameter of STAR). You do not have to set this
#'   parameter if you used the default values from STAR. Default 21.
#' @param cores Integer scalar. Number of cores to use. Default 1.
#' @param tile_width Integer scalar. The genome will be partitioned into tiles
#'   of this size. The reads in the `bam` file within each tile will be consecutively
#'   imported (or in parallel if `cores` > 1). Default 1e7.
#'
#'@return Data.frame with the coordinates of the identified novel exons. Each
#'  row in the data.frame is a predicted novel exon. The columns are the
#'  chromosome (`seqnames`), the end of the upstream exon (`lend`) in the
#'  transcript, the start and end of the novel exon (`start` and `end`), the
#'  start of the downstream exon (`rstart`) in the transcript and the strand
#'  (`strand`). The next three columns are the number of reads supporting each
#'  of the two splice junctions that define the novel exon: `unique_left` is the
#'  number of reads supporting the SJ from `lend` to `start` and `unique_right`
#'  is the number of supporting reads for the SJ from `end` to `rstart`.
#'  `min_reads` is the minimum of the two. The last column is the ID (`ID`) of
#'  the novel exon, i.e. a number between 1 and the total number of predictions.
#'
#'@importFrom data.table fread
#'@import GenomicRanges
#'@import IRanges
#'@importFrom S4Vectors queryHits subjectHits
#'@importFrom dplyr full_join mutate filter distinct
#'@importFrom GenomeInfoDb seqlevels seqlevelsInUse
#'@importFrom parallel mclapply
#'@importFrom magrittr "%>%"
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
#' ## Stranded paired-end reads where the first read comes from the reverse
#' ## strand, e.g., Illumina TruSeq stranded mRNA protocol. This is the default.
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, lib_type = "PE", stranded = "reverse")
#' 
#' ## Unstranded single-end reads
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, lib_type = "SE", stranded = "unstranded")
#'                  
#' ## If STAR was run with parameter --outSJfilterOverhangMin (e.g.
#' ## --outSJfilterOverhangMin 30 6 6 6), specify the minimal junction overhang
#' ## with parameter overhang_min: 
#' find_novel_exons(sj_filename = sj, annotation = anno, min_unique = 1,
#'                  bam = bam, overhang_min = 6)
#'                  
find_novel_exons <- function(sj_filename, annotation, min_unique = 1,
                             gzipped = FALSE, verbose = TRUE, bam,
                             single_sj = TRUE, read_based = TRUE, 
                             yield_size = 200000, 
                             lib_type = "PE", stranded = "reverse", 
                             read_length = 101, overhang_min = 12, 
                             min_intron_size = 21, cores = 1, 
                             tile_width = 1e7) {
  
  if (!lib_type %in% c("SE", "PE")) {
    stop('Parameter lib_type has to be either "SE" or "PE"')
  }
  if (!stranded %in% c("unstranded", "forward", "reverse")) {
    stop('Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  }
  if (any(read_length <= 0, overhang_min <= 0, min_intron_size <= 0)) {
    stop('Parameters "read_length", "overhang_min" and "min_intron_size" must be >= 0.')
  }
  
  exons <- annotation[["exons"]]
  introns <- annotation[["introns"]]
  
  if (verbose) message("Reading splice junctions from file")
  
  if (gzipped) {
    sj <- fread(cmd = paste0("zcat ", sj_filename))
  } else {
    sj <- fread(sj_filename)
  }
  colnames(sj) <- c("seqnames", "start", "end", "strand", "motif", "annotated",
                    "unique", "multimapping", "maxoverhang")
  
  #strand: (0: undefined, 1: +, 2: -)
  sj$strand <- c("*", "+", "-")[sj$strand + 1]
  
  ## ============= junction prefiltering ============= #
  if (verbose) message("Prefiltering splice junctions")
  
  sj_gr <- GRanges(sj)
  
  ## unify the seqnames of the novel SJs and the annotation
  ## SJs that are not present in the annotation will be removed
  sj_gr <- sj_gr[seqnames(sj_gr) %in% 
                   c(seqnames(annotation[["exons"]]), 
                     seqnames(annotation[["introns"]]))]
  GenomeInfoDb::seqlevels(sj_gr) <- seqlevelsInUse(sj_gr)

  ## filter out all annotated junctions
  sj_ann <- subsetByOverlaps(sj_gr, introns, type = "equal")
  sj_unann <- sj_gr[!(sj_gr %in% sj_ann)]
  
  ## filter out all junctions with less than min_unique reads
  sj_unann <- sj_unann[mcols(sj_unann)$unique >= min_unique]
  
  ## filter out all novel junctions that do not touch an annotated exon on the
  ## same strand
  sj_unann <- sj_unann[(start(sj_unann) - 1) %in% end(exons) |
                         (end(sj_unann) + 1) %in% start(exons), ]
  
  ## We overlap the start and end of the SJ with the genes, if both overlap with
  ## different genes, then we remove the SJ. It possibly comes from a wrongly
  ## mapped read.
  end_sj <- sj_unann
  end(end_sj) <- end(end_sj) + 1
  start(end_sj) <- end(end_sj)
  start_sj <- sj_unann
  start(start_sj) <- start(start_sj) -1
  end(start_sj) <- start(start_sj)
  
  start_genes <- findOverlaps(start_sj, genes(annotation[["txdb"]]))
  end_genes <- findOverlaps(end_sj, genes(annotation[["txdb"]]))
  sj_genes <- findOverlaps(sj_unann, genes(annotation[["txdb"]]))
  
  ## We compare the genes at the start or end of the SJ and all the genes that
  ## are located inside the SJ
  
  ## List with overlapping genes per location
  start_genes <- split(subjectHits(start_genes), queryHits(start_genes))
  end_genes <- split(subjectHits(end_genes), queryHits(end_genes))
  sj_genes <- split(subjectHits(sj_genes), queryHits(sj_genes))
  
  ind <- sapply(seq_along(sj_unann), function(x) { 
    x <- as.character(x)
    filter_trans_splice_junctions(start_genes[[x]], end_genes[[x]], 
                                  sj_genes[[x]])
    })
  
  sj_unann <- sj_unann[ind]  
  
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
      message("Importing novel SJ reads, this might take some time...")
    if (missing(bam)) {
      stop("Please specify a BAM file with parameter bam.")
    }
    
    ## annotate which end of a splice junction is touching an annotated exon
    ## We only use the exon with seqnames that are in the list of novel SJs
    touching <- mclapply(seq_along(sj_unann), function(i) 
      sj_touching_exon(as.character(seqnames(sj_unann[i])), 
                       start(sj_unann[i]), end(sj_unann[i]), 
                       strand(sj_unann[i]), 
                       exons = exons[seqnames(exons) %in% seqnames(sj_unann)]), 
      mc.cores = cores)
    touching <- unlist(touching)
    
    sj_unann$touching <- touching
    ## remove the ambiguous sj that touch different genes on both ends
    sj_unann <- sj_unann[!is.na(touching), ]
    
    ## TODO: is this really necessary? This is faster when using 15 cores...
    if(cores > 5){  
      reads <- import_novel_sj_reads_parallel(bam, yield_size = yield_size, 
                                              sj_unann, cores = cores)
      reads <- do.call("c", reads) 
    } else {
      reads <- import_novel_sj_reads(bam, sj_unann)
    }

    sj_unann <- data.frame(sj_unann, stringsAsFactors = TRUE)
    novel_exons <- rbind(novel_exons, 
                         identify_exon_from_sj(sj_unann, reads = reads,
                                               txdb = annotation[["txdb"]], 
                                               cores = cores))
    ## remove unneeded objects
    rm(reads, sj_unann, touching)
  }
  
  if (read_based) {
    ## ============ Predict novel exons from reads with 2 junctions ==============
    if (verbose) message("Predict novel exons from reads with 2 junctions")
    
    if (missing(bam)) {
      stop("Please specify a BAM file with parameter bam.")
    }

    junc_reads <- filter_junction_reads(bam, lib_type = lib_type, 
                                        stranded = stranded, cores = cores, 
                                        yield_size = yield_size, 
                                        tile_width = tile_width)
    ## unlist the list of GAlignments objects
    junc_reads <- do.call(c, junc_reads)  

    if (length(junc_reads) == 0){
      stop(paste0("Reading the BAM file", bam, 
                  " resulted in an empty GAlignments object."))
    }
    
    read_pred <- predict_jr_exon(junc_reads, annotation)
    
    if (exists("novel_exons", inherits = FALSE)) {
      ## make sure the factors have the same levels
      novel_exons <- merge(novel_exons, read_pred, all = TRUE)
      novel_exons$seqnames <- droplevels(novel_exons$seqnames)
    } else {
      novel_exons <- read_pred
    }
    
    ## ======= Predict novel exons from read pairs with each one junction =======
    if (lib_type == "PE") {
      if (verbose) message("Predict novel exons from read pairs with each 1 junction")
      read_pair_pred <- predict_jrp_exon(junc_reads, annotation, 
                                         read_length = read_length, 
                                         overhang_min = overhang_min, 
                                         min_intron_size = min_intron_size)
      novel_exons <- merge(novel_exons, read_pair_pred, all = TRUE)
      novel_exons$seqnames <- droplevels(novel_exons$seqnames)
    }
  }
  
  ## Filter the table with novel exons and remove all duplicate predictions
  novel_exons <- distinct(novel_exons)
  
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
  
  ## We remove all predicted cassette exons where one of the junctions is not in
  ## SJ.out.tab, because these exons are most likely wrong (this only happens if
  ## we predict exons from reads).
  novel_exons <- novel_exons %>% 
    filter(!(!is.na(lend) & is.na(unique_left))) %>%
    filter(!(!is.na(rstart) & is.na(unique_right)))
  
  ## take the minimum read coverage of both junctions
  novel_exons$min_reads <- pmin(novel_exons$unique_left,
                                novel_exons$unique_right, na.rm = TRUE)
  novel_exons$ID <- 1:nrow(novel_exons)
  
  ## Filter all predictions with less than `min_unique` reads supporting both
  ## SJs
  novel_exons <- novel_exons %>% filter(min_reads >= min_unique)
  novel_exons
}
