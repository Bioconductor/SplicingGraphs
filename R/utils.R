### =========================================================================
### Some low-level internal utilities
### -------------------------------------------------------------------------


### Raises an error if the ranges in GRanges object 'x' are not all on the
### same reference sequence and strand.
### Returns "+", "-", or "*" (returns "*" if 'x' is empty).
commonStrand.GRanges <- function(x, what="ranges in 'x'")
{
    if (!is(x, "GRanges"))
        stop("'x' must be a GRanges object")
    ## An arbitrary choice. Shouldn't have any impact on what happens later on.
    if (length(x) == 0L)
        return("*")
    x_strand <- strand(x)
    x_strand_runValue <- runValue(x_strand)
    if (nrun(seqnames(x)) != 1L || length(x_strand_runValue) != 1L)
        stop("all the ", what, " must be on the same ",
             "reference sequence and strand")
    as.vector(x_strand_runValue[1L])
}

### Returns a factor, NOT a factor-Rle!
commonSeqnames.RleList <- function(x, errmsg="internal error")
{
    if (!is(x, "RleList"))
        stop("'x' must be a RleList object")
    x_runValue <- runValue(x)
    x_runValue_eltlens <- elementLengths(x_runValue)
    if (!all(x_runValue_eltlens <= 1L))
        stop(errmsg)
    idx0 <- which(x_runValue_eltlens == 0L)
    idx0_len <- length(idx0)
    if (idx0_len == 0L)
        return(x_runValue@unlistData)
    no_val <- factor(NA_character_, levels=levels(x_runValue@unlistData))
    x_runValue[idx0] <- as(rep.int(no_val, idx0_len), "List")
    IRanges:::decodeRle(x_runValue@unlistData)
}

### Returns a factor, NOT a factor-Rle!
commonStrand.RleList <- function(x, errmsg="internal error")
{
    if (!is(x, "RleList"))
        stop("'x' must be a RleList object")
    x_runValue <- runValue(x)
    x_runValue_eltlens <- elementLengths(x_runValue)
    if (!all(x_runValue_eltlens <= 1L))
        stop(errmsg)
    idx0 <- which(x_runValue_eltlens == 0L)
    idx0_len <- length(idx0)
    if (idx0_len == 0L)
        return(x_runValue@unlistData)
    no_val <- strand("*")
    x_runValue[idx0] <- as(rep.int(no_val, idx0_len), "List")
    ans <- x_runValue@unlistData
    factor(levels(strand())[ans], levels=levels(strand()))
}

commonSeqnames.GRangesList <- function(x)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    errmsg <- c("some top-level elements in 'x' contain ranges that ",
                "are not all on the same reference sequence and strand")
    Rle(commonSeqnames.RleList(seqnames(x), errmsg=errmsg))
}

commonStrand.GRangesList <- function(x)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    x_seqnames <- seqnames(x)
    x_seqnames_runValue <- runValue(x_seqnames)
    x_seqnames_runValue_eltlens <- elementLengths(x_seqnames_runValue)
    errmsg <- c("some top-level elements in 'x' contain ranges that ",
                "are not all on the same reference sequence and strand")
    if (!all(x_seqnames_runValue_eltlens <= 1L))
        stop(errmsg)
    Rle(commonStrand.RleList(strand(x), errmsg=errmsg))
}

