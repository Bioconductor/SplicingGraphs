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

commonStrand.GRangesList <- function(x)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    x_seqnames <- seqnames(x)
    x_seqnames_runValue <- runValue(x_seqnames)
    x_seqnames_runValue_eltlens <- elementLengths(x_seqnames_runValue)
    x_strand <- strand(x)
    x_strand_runValue <- runValue(x_strand)
    x_strand_runValue_eltlens <- elementLengths(x_strand_runValue)
    if (!all(x_seqnames_runValue_eltlens <= 1L) ||
        !all(x_strand_runValue_eltlens <= 1L))
        stop("some top-level elements in 'x' contain ranges that ",
             "are not all on the same reference sequence and strand")
    idx <- which(x_strand_runValue_eltlens == 0L)
    idx_len <- length(idx)
    if (idx_len == 0L)
        return(Rle(x_strand_runValue@unlistData))
    x_strand_runValue[idx] <- as(rep.int(strand("*"), idx_len), "List")
    ans <- x_strand_runValue@unlistData
    Rle(factor(levels(strand())[ans], levels=levels(strand())))
}

