###

load_TSPC_gene_model <- function(models_path, check.transcripts=TRUE)
{
    models <- read.table(models_path, stringsAsFactors=FALSE)
    stopifnot(ncol(models) == 3L)  # sanity check
    tmp1 <- strsplit(models[[1L]], ":", fixed=TRUE)
    stopifnot(all(elementLengths(tmp1) == 2L))  # sanity check
    tmp1 <- unlist(tmp1, use.name=FALSE)
    exons_seqnames <- tmp1[c(TRUE, FALSE)]
    exons_ranges <- tmp1[c(FALSE, TRUE)]
    tmp2 <- strsplit(exons_ranges, "-", fixed=TRUE)
    stopifnot(all(elementLengths(tmp2) == 2L))  # sanity check
    tmp2 <- as.integer(unlist(tmp2, use.name=FALSE))
    stopifnot(!any(is.na(tmp2)))  # sanity check
    exons_start <- tmp2[c(TRUE, FALSE)]
    exons_end <- tmp2[c(FALSE, TRUE)]
    exons_ranges <- IRanges(exons_start, exons_end)
    exon_id <- rep.int(NA_integer_, length(exons_ranges))
    exon_rank <- IRanges:::fancy_mseq(runLength(Rle(models[[2L]])))
    unlisted_ans <- GRanges(seqnames=exons_seqnames,
                            ranges=exons_ranges,
                            exon_id=exon_id,
                            exon_name=models[[3L]],
                            exon_rank=exon_rank)
    ans <- split(unlisted_ans, models[[2L]])
    if (check.transcripts)
        stopifnot(all(isNormal(ranges(ans))))
    strand(ans@unlistData) <- "+"
    ans
}

get_TSPC_models_path <- function(subdir_path)
{
    if (!isSingleString(subdir_path))
        stop("'subdir_path' must be a single string")
    SUFFIX <- "_models.txt"
    filenames <- list.files(subdir_path)
    stop <- nchar(filenames)
    start <- stop - nchar(SUFFIX) + 1L
    suffixes <- substr(filenames, start, stop)
    models_filename <- filenames[suffixes == SUFFIX]
    if (length(models_filename) != 1L)
        stop("found more than one models file in ", subdir_path)
    file.path(subdir_path, models_filename)
}

make_TSPC_SplicinGraphs <- function(subdir_paths)
{
    gene_list <- lapply(subdir_paths,
        function(subdir_path) {
            models_path <- get_TSPC_models_path(subdir_path)
            message("Reading ", models_path, " ... ", appendLF=FALSE)
            gene <- load_TSPC_gene_model(models_path)
            message("OK")
            gene
        })
    suppressWarnings(genes <- do.call(c, unname(gene_list)))
    genes_seqlevels <- seqlevels(genes)
    seqlevels(genes) <- genes_seqlevels[order(makeSeqnameIds(genes_seqlevels))]
    grouping <- rep.int(basename(subdir_paths), elementLengths(gene_list))
    SplicingGraphs(genes, grouping=grouping, min.ntx=1L)
}

get_TSPC_sample_names <- function(subdir_paths)
{
    SUFFIX <- ".bam"
    sample_names <- lapply(subdir_paths,
        function(subdir_path) {
            filenames <- list.files(subdir_path)
            stop <- nchar(filenames)
            start <- stop - nchar(SUFFIX) + 1L
            suffixes <- substr(filenames, start, stop)
            bam_filenames <- filenames[suffixes == SUFFIX]
            prefix <- paste0(basename(subdir_path), "-")
            start <- nchar(prefix) + 1L
            stop <- nchar(bam_filenames) - nchar(SUFFIX)
            ans <- substr(bam_filenames, start, stop)
            stopifnot(!anyDuplicated(ans))
            ans
        })
    unique(unlist(sample_names, use.names=FALSE))
}

get_TSPC_bam_path <- function(subdir_path, sample_name)
{
    if (!isSingleString(subdir_path))
        stop("'subdir_path' must be a single string")
    if (!isSingleString(sample_name))
        stop("'sample_name' must be a single string")
    SUFFIX <- ".bam"
    prefix <- paste0(basename(subdir_path), "-")
    bam_filename <- paste0(prefix, sample_name, SUFFIX)
    file.path(subdir_path, bam_filename)
}

### 7 genes with samples: BAI1, CYB561, DAPL1, ITGB8, LGSN, MKRN3, ST14.
### 54 samples:
###   - 34 samples are single-end for all genes.
###   - 12 samples are paired-end for all genes.
###   - 8 samples are single-end for some genes and paired-end for the others.
###     Those samples are: SOC_11199_480891, SOC_2774_227362, SOC_5973_278042,
###     SOC_7244_348592, SOC_7637_361542, SOC_7777_371281, SOC_8904_437502,
###     and SOC_9547_467919.
is_paired_end_TSPC_sample <- function(sample_names, subdir_path)
{
    if (!isSingleString(subdir_path))
        stop("'subdir_path' must be a single string")
    library(Rsamtools)
    flag0 <- scanBamFlag(#isProperPair=TRUE,
                         isNotPrimaryRead=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0, what="flag")
    sapply(sample_names,
        function(sample_name) {
            bam_filepath <- get_TSPC_bam_path(subdir_path, sample_name)
            res <- scanBam(bam_filepath, use.names=TRUE, param=param0)
            flag <- bamFlagTest(res[[1L]]$flag, "isPaired")
            if (all(flag))
                return(TRUE)
            if (any(flag))
                return(NA)
            FALSE
        })
}

### Returns the reads as a GRangesList object.
load_TSPC_sample_reads <- function(sample_name, subdir_paths)
{
    library(Rsamtools)
    flag0 <- scanBamFlag(#isProperPair=TRUE,
                         isNotPrimaryRead=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0, what=c("flag", "mapq"))
    reads_list <- lapply(subdir_paths,
        function(subdir_path) {
            bam_filepath <- get_TSPC_bam_path(subdir_path, sample_name)
            if (!file.exists(bam_filepath)) {
                return(GRangesList())
            }
            gal <- readGappedAlignments(bam_filepath, use.names=TRUE,
                                        param=param0)
            is_paired <- bamFlagTest(mcols(gal)$flag, "isPaired")
            if (!any(is_paired)) {
                ## The aligner reported 2 *primary* alignments for single-end
                ## read s100208_3_83_5646_14773 in file
                ## BAI1-SOC_5991_294171.bam, which doesn't make sense. However,
                ## the reported mapping quality for those 2 alignments is 3
                ## which is very low. So let's get rid of alignments that have
                ## a quality <= 3.
                mapq <- mcols(gal)$mapq
                gal <- gal[is.na(mapq) | mapq > 3]
                ans <- grglist(gal, order.as.in.query=TRUE)
            } else {
                stopifnot(all(is_paired))
                galp <- readGappedAlignmentPairs(bam_filepath, use.names=TRUE,
                                                 param=param0)
                ans <- grglist(galp, order.as.in.query=TRUE)
           }
           ans
        })
    ans <- do.call(c, unname(reads_list))
    stopifnot(!anyDuplicated(names(ans)))
    ans
}

assign_TSPC_reads <- function(sg, subdir_paths)
{
    sample_names <- get_TSPC_sample_names(subdir_paths)
    nsample <- length(sample_names)
    for (i in seq_len(nsample)) {
        sample_name <- sample_names[[i]]
        message("Assign reads from sample ", sample_name,
                " (", i, "/", nsample, ") ... ", appendLF=FALSE)
        reads <- load_TSPC_sample_reads(sample_name, subdir_paths)
        sg <- assignReads(sg, reads, sample.name=sample_name)
        message("OK")
    }
    sg
}

