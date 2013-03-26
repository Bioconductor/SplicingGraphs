### =========================================================================
### Some low-level internal utilities
### -------------------------------------------------------------------------


### Raises an error if the ranges in GRanges object 'x' are not all on the
### same reference sequence and strand.
### Returns "+", "-", or "*" (returns "+" if 'x' is empty).
commonStrand <- function(x, what="ranges in 'x'")
{
    if (!is(x, "GRanges"))
        stop("'x' must be a GRanges object")
    ## An arbitrary choice. Shouldn't have any impact on what happens later on.
    if (length(x) == 0L)
        return("+")
    x_strand <- strand(x)
    x_strand_runValue <- runValue(x_strand)
    if (nrun(seqnames(x)) != 1L || length(x_strand_runValue) != 1L)
        stop("all the ", what, " must be on the same ",
             "reference sequence and strand")
    as.vector(x_strand_runValue[1L])
}

