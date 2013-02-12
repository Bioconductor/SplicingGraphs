###

loadModels <- function(models_path, check.transcripts=TRUE)
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
    unlisted_ans <- GRanges(seqnames=exons_seqnames,
                            ranges=IRanges(exons_start, exons_end),
                            exon_name=models[[3L]])
    ans <- split(unlisted_ans, models[[2L]])
    if (check.transcripts)
        stopifnot(all(isNormal(ranges(ans))))
    strand(ans@unlistData) <- "+"
    ans
}

### It's questionable whether this does the right thing on paired-end reads.
### I guess not...
makeSgdfWithHits <- function(grl, ex_by_tx2)
{
    ov0 <- findOverlaps(grl, ex_by_tx2, ignore.strand=TRUE)
    ovenc0 <- encodeOverlaps(grl, ex_by_tx2, hits=ov0,
                             flip.query.if.wrong.strand=TRUE)
    ov0_is_comp <- isCompatibleWithSplicing(ovenc0)
    ov1 <- ov0[ov0_is_comp]
    ex_by_tx2 <- assignSubfeatureHits(grl, ex_by_tx2, ov1, ignore.strand=TRUE)
    in_by_tx2 <- psetdiff(range(ex_by_tx2), ex_by_tx2)
    in_by_tx2 <- assignSubfeatureHits(grl, in_by_tx2, ov1, ignore.strand=TRUE)
    Sgdf(ex_by_tx2, inbytx=in_by_tx2)
}

makeTSPCsgdf <- function(subdir_path)
{
    subdir_basename <- basename(subdir_path)
    filenames <- list.files(subdir_path)
    filenames_nchar <- nchar(filenames)

    ## Load the gene model.
    suffixes <- substr(filenames, filenames_nchar-10L, filenames_nchar)
    models_filename <- filenames[suffixes == "_models.txt"]
    models_path <- file.path(subdir_path, models_filename)
    message("Reading ", models_path, " ... ", appendLF=FALSE)
    ex_by_tx <- loadModels(models_path)
    message("OK")

    ## Compute the splicing graph.
    ex_by_tx2 <- splicingGraphs(ex_by_tx)
    ans <- Sgdf(ex_by_tx2)

    ## Find the BAM files.
    suffixes <- substr(filenames, filenames_nchar-3L, filenames_nchar)
    bam_filenames <- filenames[suffixes == ".bam"]
    nbam <- length(bam_filenames)
    if (nbam == 0L)
        return(ans)

    ## Load and process the BAM files.
    prefixes <- substr(bam_filenames, 1L, nchar(subdir_basename)+1L)
    stopifnot(all(prefixes == paste0(subdir_basename, "-")))
    sample_names <- substr(bam_filenames, nchar(prefixes)+1L,
                                          nchar(bam_filenames)-4L)
    flag0 <- scanBamFlag(#isProperPair=TRUE,
                         isNotPrimaryRead=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0, what=c("flag", "mapq"))
    X <- seq_len(nbam)
    names(X) <- sample_names
    nhits <- sapply(X, function(i) {
        bam_filename <- bam_filenames[i]
        sample_name <- sample_names[i]
        message("Processing BAM file ", i, "/", nbam,
                " (sample ", sample_name, ") ... ", appendLF=FALSE)
        bam_filepath <- file.path(subdir_path, bam_filename)
        gal <- readGappedAlignments(bam_filepath, use.names=TRUE,
                                    param=param0)
        is_paired <- bamFlagTest(mcols(gal)$flag, "isPaired")
        if (!any(is_paired)) {
            ## The aligner reported 2 *primary* alignments for single-end read
            ## s100208_3_83_5646_14773 in file BAI1-SOC_5991_294171.bam, which
            ## doesn't make sense. However, the reported mapping quality for
            ## those 2 alignments is 3 which is very low. So let's get rid of
            ## alignments that have a quality <= 3.
            mapq <- mcols(gal)$mapq
            gal <- gal[is.na(mapq) | mapq > 3]
            grl <- grglist(gal, order.as.in.query=TRUE)
        } else {
            stopifnot(all(is_paired))
            galp <- readGappedAlignmentPairs(bam_filepath, use.names=TRUE,
                                             param=param0)
            grl <- grglist(galp, order.as.in.query=TRUE)
        }
        sgdf <- makeSgdfWithHits(grl, ex_by_tx2)
        message("OK")
        sgdf[ , "nhits"]
    })
    cbind(ans, DataFrame(nhits))
}

