### =========================================================================
### Some low-level internal utilities
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the common sequence name or strand from various kinds of objects
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Predefined edge metadata columns
###

### Edge metadata columns that are considered to be exon attributes (note
### that we include the "start_SSid" and "end_SSid" cols). Those columns are
### the 5 first inner metadata columns of the GRangesList object containing
### the exons grouped by transcripts returned by unlist() when called on a
### SplicingGraphs object.
EXON_MCOLS <- c("exon_id", "exon_name", "exon_rank", "start_SSid", "end_SSid")

### All edge metadata columns.
ALL_EDGE_MCOLS <- c("sgedge_id", "from", "to", "ex_or_in", "tx_id", EXON_MCOLS)

### Subset of 'ALL_EDGE_MCOLS' made of those columns that are considered
### invariant i.e. the values in them associated with the same sgedge_id
### (global edge id) should be the same. Note that we also include the
### "sgedge_id" col itself.
INVARIANT_EDGE_MCOLS <- c("sgedge_id", "from", "to", "ex_or_in",
                          "start_SSid", "end_SSid")

EX_OR_IN_LEVELS2 <- c("ex", "in", "", "mixed")
EX_OR_IN_LEVELS <- EX_OR_IN_LEVELS2[-4L]

valid_exon_mcolnames <- function(colnames)
{
    nb_exon_mcols <- length(EXON_MCOLS)
    if (identical(head(colnames, n=nb_exon_mcols), EXON_MCOLS))
        return(NULL)
    msg <- c("the first ", nb_exon_mcols, " exon-level metadata columns ",
             "must be: ", paste0(EXON_MCOLS, collapse=", "))
    paste0(msg, collapse="")
}

check_exon_mcolnames <- function(colnames)
{
    msg <- valid_exon_mcolnames(colnames)
    if (!is.null(msg))
        stop(msg)
}

check_all_edge_mcolnames <- function(colnames)
{
    stopifnot(identical(head(colnames, n=length(ALL_EDGE_MCOLS)),
                        ALL_EDGE_MCOLS))
}

get_index_of_group_of_mcols <- function(colnames,
                                        with.exon.mcols, with.hits.mcols)
{
    ans <- integer(0)
    if (!with.exon.mcols) {
        idx <- match(EXON_MCOLS, colnames)
        ans <- c(ans, idx)
    }
    if (!with.hits.mcols) {
        idx <- grep("hits$", colnames)
        ans <- c(ans, idx)
    }
    ans
}

get_index_of_invariant_edge_mcols <- function(colnames)
{
    idx <- match(INVARIANT_EDGE_MCOLS, colnames)
    idx[!is.na(idx)]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### unlistAndSplit()
###
### Example:
###
###   > x <- SimpleList(A=4:5, B=letters[1:4], C=NULL, D=1:2, E=-2:0, F=TRUE)
###   > f <- c("y", "x", "x", "x", "y", "x")
###   > unlistAndSplit(x, f)
###   CharacterList of length 2
###   [["x"]] a b c d 1 2 TRUE
###   [["y"]] 4 5 -2 -1 0
###   > unlistAndSplit(x, f)[[1]]
###        B      B      B      B      D      D      F 
###      "a"    "b"    "c"    "d"    "1"    "2" "TRUE" 
###
### Should work on any vector-like object and act as an endomorphism on a
### CompressedList object. On an atomic vector (on which 'unlist()' is a
### no-op), should be equivalent to 'splitAsList(x, f)'.
### TODO: Maybe move this to IRanges and expose to the user.
unlistAndSplit <- function(x, f, drop=FALSE)
{
    if (length(f) != length(x))
        stop("'x' and 'f' must have the same length")
    x2 <- unlist(x)
    f2 <- rep.int(f, elementLengths(x))
    splitAsList(x2, f2, drop=drop)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### uninformative_sgnodes()
###

### TODO: Rename this hasDuplicates() and move it to IRanges.
.has_duplicates <- function(x) (duplicated(x) | duplicated(x, fromLast=TRUE))

.singletons <- function(x) x[!.has_duplicates(x)]

uninformative_sgnodes <- function(from, to)
{
    from1_nodes <- .singletons(from)
    to1_nodes <- .singletons(to)
    intersect(from1_nodes, to1_nodes)
}

