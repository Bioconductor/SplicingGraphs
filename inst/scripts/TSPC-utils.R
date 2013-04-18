### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load the gene models
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
    ## The '_models.txt' files use 0-based starts!
    exons_start <- tmp2[c(TRUE, FALSE)] + 1L
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Make the TSPC SplicinGraphs object
###

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
    suppressWarnings(ex_by_tx <- do.call(c, unname(gene_list)))
    ex_by_tx_seqlevels <- seqlevels(ex_by_tx)
    seq_rank <- makeSeqnameIds(ex_by_tx_seqlevels)
    seqlevels(ex_by_tx) <- ex_by_tx_seqlevels[order(seq_rank)]
    grouping <- rep.int(basename(subdir_paths), elementLengths(gene_list))
    SplicingGraphs(ex_by_tx, grouping=grouping, min.ntx=1L)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load the sample reads
###

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

### BAM status:
###   ".": BAM file doesn't exist;
###   "0": file is empty (no alignments);
###   "s": single-end;
###   "p": paired-end;
###   "m": mixed single-/paired-end.
get_TSPC_bam_status <- function(subdir_path, sample_name)
{
    bam_filepath <- get_TSPC_bam_path(subdir_path, sample_name)
    if (!file.exists(bam_filepath))
        return(".")
    library(Rsamtools)
    flag0 <- scanBamFlag(#isProperPair=TRUE,
                         isNotPrimaryRead=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0, what="flag")
    res <- scanBam(bam_filepath, use.names=TRUE, param=param0)
    stopifnot(length(res) == 1L)
    flag <- res[[1L]]$flag
    nb_rec <- length(flag)
    if (nb_rec == 0L)
        return("0")
    nb_paired <- sum(bamFlagTest(flag, "isPaired"))
    if (nb_paired == 0L)
        return("s")
    if (nb_paired == nb_rec)
        return("p")
    "m"
}

### Returns a matrix with 1 row per path in 'subdir_paths', and 1 col per
### sample in 'sample_names'.
make_TSPC_bam_status_matrix <- function(subdir_paths, sample_names)
{
    if (!is.character(subdir_paths))
        stop("'subdir_paths' must be a character vector")
    if (!is.character(sample_names))
        stop("'sample_names' must be a character vector")
    ans <- sapply(sample_names,
        function(sample_name)
            sapply(subdir_paths, get_TSPC_bam_status, sample_name))
    rownames(ans) <- basename(subdir_paths)
    ans
}

get_TSPC_bam_gaprate <- function(subdir_path, sample_name)
{
    bam_filepath <- get_TSPC_bam_path(subdir_path, sample_name)
    if (!file.exists(bam_filepath))
        return(NA_real_)
    library(Rsamtools)
    flag0 <- scanBamFlag(#isProperPair=TRUE,
                         isNotPrimaryRead=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0, what="cigar")
    res <- scanBam(bam_filepath, use.names=TRUE, param=param0)
    stopifnot(length(res) == 1L)
    cigar <- res[[1L]]$cigar
    nb_rec <- length(cigar)
    length(grep("N", cigar, fixed=TRUE)) / nb_rec
}

### Returns a matrix with 1 row per path in 'subdir_paths', and 1 col per
### sample in 'sample_names'.
make_TSPC_bam_gaprate_matrix <- function(subdir_paths, sample_names)
{
    if (!is.character(subdir_paths))
        stop("'subdir_paths' must be a character vector")
    if (!is.character(sample_names))
        stop("'sample_names' must be a character vector")
    ans <- sapply(sample_names,
        function(sample_name)
            sapply(subdir_paths, get_TSPC_bam_gaprate, sample_name))
    rownames(ans) <- basename(subdir_paths)
    ans
}

### Returns a GAlignments or GAlignmentPairs object, or NULL if the file
### doesn't exist or is empty (no alignments).
load_TSPC_bam <- function(subdir_path, sample_name, gapped.reads.only=FALSE)
{
    if (!isTRUEorFALSE(gapped.reads.only))
        stop("'gapped.reads.only' must be TRUE or FALSE")
    bam_status <- get_TSPC_bam_status(subdir_path, sample_name)
    message(bam_status, appendLF=FALSE)
    if (bam_status %in% c(".", "0"))
        return(NULL)
    bam_filepath <- get_TSPC_bam_path(subdir_path, sample_name)
    if (bam_status == "m")
        stop("file ", bam_filepath, " contains a mix of single- ",
             "and paired-end reads")
    library(Rsamtools)
    flag0 <- scanBamFlag(#isProperPair=TRUE,
                         isNotPrimaryRead=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0, what="mapq")
    if (bam_status == "p") {
        reads <- readGAlignmentPairs(bam_filepath, use.names=TRUE,
                                     param=param0)
        if (!gapped.reads.only)
            return(reads)
        keep_idx <- which(grepl("N", cigar(first(reads)), fixed=TRUE) |
                          grepl("N", cigar(last(reads)), fixed=TRUE))
        return(reads[keep_idx])
    }
    reads <- readGAlignments(bam_filepath, use.names=TRUE, param=param0)
    ## The aligner reported 2 *primary* alignments for single-end
    ## read s100208_3_83_5646_14773 in file BAI1-SOC_5991_294171.bam, which
    ## doesn't make sense. However, the reported mapping quality for those
    ## 2 alignments is 3 which is very low. So let's get rid of alignments
    ## that have a quality <= 3.
    mapq <- mcols(reads)$mapq
    lowmapq_idx <- which(!is.na(mapq) & mapq <= 3)
    if (length(lowmapq_idx) != 0L) {
        message("|lowmapq:", length(lowmapq_idx), appendLF=FALSE)
        reads <- reads[-lowmapq_idx]
    }
    if (!gapped.reads.only)
        return(reads)
    keep_idx <- which(grepl("N", cigar(reads), fixed=TRUE))
    reads[keep_idx]
}

### Returns a GAlignments or GAlignmentPairs object.
load_TSPC_sample_reads <- function(subdir_paths, sample_name,
                                   gapped.reads.only=FALSE)
{
    reads_list <- lapply(subdir_paths,
        function(subdir_path) {
            message("<", basename(subdir_path), "|", appendLF=FALSE)
            reads <- load_TSPC_bam(subdir_path, sample_name,
                                   gapped.reads.only=gapped.reads.only)
            if (!is.null(reads))
                message("|", length(reads), appendLF=FALSE)
            message("> ", appendLF=FALSE)
            reads
        })
    empty_idx <- which(elementLengths(reads_list) == 0L)
    if (length(empty_idx) != 0L)
        reads_list <- reads_list[-empty_idx]
    reads <- do.call(c, unname(reads_list))
    stopifnot(!anyDuplicated(names(reads)))
    reads
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compute matrix of "Gapped Read Compatibility Ratio"
###

get_TSPC_grcr <- function(sg, subdir_path, sample_name)
{
    suppressMessages(reads <- load_TSPC_bam(subdir_path, sample_name,
                                            gapped.reads.only=TRUE))
    if (is.null(reads))
        return(NA_real_)
    ex_by_tx <- sg[[basename(subdir_path)]]
    ov <- findCompatibleOverlaps(reads, ex_by_tx)
    sum(countQueryHits(ov) != 0L) / length(reads)
}

make_TSPC_grcr_matrix <- function(sg, subdir_paths, sample_names)
{
    if (!is(sg, "SplicingGraphs"))
        stop("'sg' must be a SplicingGraphs object")
    if (!is.character(subdir_paths))
        stop("'subdir_paths' must be a character vector")
    if (!is.character(sample_names))
        stop("'sample_names' must be a character vector")
    ans <- sapply(sample_names,
        function(sample_name)
            sapply(subdir_paths, get_TSPC_grcr, sg=sg, sample_name=sample_name))
    rownames(ans) <- basename(subdir_paths)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Assign the TSPC reads to the SplicingGraphs object
###

assign_TSPC_reads <- function(sg, subdir_paths, sample_names,
                              gapped.reads.only=FALSE)
{
    if (!is(sg, "SplicingGraphs"))
        stop("'sg' must be a SplicingGraphs object")
    if (!is.character(subdir_paths))
        stop("'subdir_paths' must be a character vector")
    if (!is.character(sample_names))
        stop("'sample_names' must be a character vector")
    nsample <- length(sample_names)
    for (i in seq_len(nsample)) {
        sample_name <- sample_names[[i]]
        message("Assign reads from sample ", sample_name,
                " (", i, "/", nsample, ") ... ", appendLF=FALSE)
        reads <- load_TSPC_sample_reads(subdir_paths, sample_name,
                                        gapped.reads.only=gapped.reads.only)
        sg <- assignReads(sg, reads, sample.name=sample_name)
        message("OK")
    }
    sg
}

