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

