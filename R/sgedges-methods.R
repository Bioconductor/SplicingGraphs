### =========================================================================
### "sgedges" (and related) methods
### -------------------------------------------------------------------------


.get_sgnodes_from_txpaths <- function(txpaths)
{
    SSids <- unique(unlist(txpaths, use.names=FALSE))
    c("R", sort(SSids), "L")
}

.get_sgnodes_from_sgedges <- function(sgedges)
{
    from <- sgedges[ , "from"]
    to <- sgedges[ , "to"]
    SSids <- as.integer(setdiff(c(from, to), c("R", "L")))
    c("R", sort(SSids), "L")
}

make_matrix_from_txpaths <- function(txpaths)
{
    sgnodes <- .get_sgnodes_from_txpaths(txpaths)
    ans_nrow <- length(txpaths)
    ans_ncol <- length(sgnodes)
    ans_dimnames <- list(names(txpaths), sgnodes)
    ans <- matrix(FALSE , nrow=ans_nrow, ncol=ans_ncol, dimnames=ans_dimnames)
    ans[ , 1L] <- ans[ , ans_ncol] <- TRUE
    i <- cbind(rep.int(seq_along(txpaths), elementLengths(txpaths)),
               unlist(txpaths, use.names=FALSE) + 1L)
    ans[i] <- TRUE
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### txpaths() accessor
###
### Gets the splicing paths.
### Returns them in a named IntegerList with 1 top-level element per
### transcript in the specified gene. Each top-level element 'txpaths[[i]]'
### contains the splicing site ids for the i-th transcript.
###

setGeneric("txpaths", signature="x",
    function(x, as.matrix=FALSE) standardGeneric("txpaths")
)

### Should return a CompressedIntegerList.
setMethod("txpaths", "SplicingGraphs",
    function(x, as.matrix=FALSE)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        if (!isTRUEorFALSE(as.matrix))
            stop("'as.matrix' must be TRUE or FALSE")
        ans <- mcols(unlist(x, use.names=FALSE))[ , "txpaths"]
        if (as.matrix)
            ans <- make_matrix_from_txpaths(ans)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UATXHcount() accessor
###

setGeneric("UATXHcount", signature="x",
    function(x) standardGeneric("UATXHcount")
)

### Should return an integer vector or a NULL.
setMethod("UATXHcount", "SplicingGraphs",
    function(x)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        mcols(unlist(x, use.names=FALSE))[["UATXHcount"]]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgedges() extractor
###
### Returns the splicing graph in a DataFrame with 1 row per edge.
###

### 'txpaths' must be an IntegerList containing all the splicing paths (1 per
### transcript) for a given gene. Should have been obtained thru the txpaths()
### accessor. Returns a 4-col (or 5-col if 'UATXHcount' is supplied) data.frame
### representing the splicing graph.
.make_sgedges0_from_txpaths <- function(txpaths, UATXHcount=NULL)
{
    if (!is.null(UATXHcount)) {
        if (!is.integer(UATXHcount))
            stop("'UATXHcount' must be an integer vector or NULL")
        if (length(UATXHcount) != length(txpaths))
            stop("when not NULL, 'UATXHcount' must have ",
                 "the same length as 'txpaths'")
    }
    sgedges0s <- lapply(seq_along(txpaths),
                        function(i) {
                            txpath <- txpaths[[i]]
                            txpath_len <- length(txpath)
                            if (txpath_len %% 2L != 0L)
                                stop("some paths in 'txpaths' contain ",
                                     "an odd number of splicing site ids")
                            from <- c("R", txpath)
                            to <- c(txpath, "L")
                            nexons <- txpath_len %/% 2L
                            if (nexons == 0L) {
                                ex_or_in <- EX_OR_IN_LEVELS[3L]
                            } else {
                                nintrons <- nexons - 1L
                                ex_or_in <- c(EX_OR_IN_LEVELS[3L],
                                              rep.int(EX_OR_IN_LEVELS[1:2],
                                                      nintrons),
                                              EX_OR_IN_LEVELS[1L],
                                              EX_OR_IN_LEVELS[3L])
                            }
                            ex_or_in <- factor(ex_or_in,
                                               levels=EX_OR_IN_LEVELS)
                            data.frame(from=from,
                                       to=to,
                                       ex_or_in=ex_or_in,
                                       stringsAsFactors=FALSE)
                        })
    nedges_per_tx <- sapply(sgedges0s, nrow)
    sgedges0 <- do.call(rbind, sgedges0s)
    tx_id <- names(txpaths)
    if (is.null(tx_id))
        tx_id <- seq_along(txpaths)
    tx_id <- rep.int(factor(tx_id, levels=tx_id), nedges_per_tx)
    sgedges0$tx_id <- tx_id
    if (!is.null(UATXHcount))
        sgedges0$UATXHcount <- rep.int(UATXHcount, nedges_per_tx)
    sgedges0
}

### Collapse the duplicated edges in 'sgedges0' into a DataFrame.
### We use a DataFrame instead of a data.frame because we want to store
### the tx_id col in a CompressedFactorList (even though this container
### doesn't formally exist and a CompressedIntegerList is actually what's
### being used).
.make_sgedges_from_sgedges0 <- function(sgedges0, ex_hits=NULL, in_hits=NULL)
{
    from <- sgedges0[ , "from"]
    to <- sgedges0[ , "to"]
    ex_or_in <- sgedges0[ , "ex_or_in"]
    tx_id <- sgedges0[ , "tx_id"]
    edges <- paste(from, to, sep="~")
    sm <- match(edges, edges)
    if (!all(ex_or_in == ex_or_in[sm]))
        stop("invalid splicing graph")
    is_not_dup <- sm == seq_along(sm)
    sgedges <- DataFrame(sgedges0[is_not_dup, , drop=FALSE])
    sgedges$tx_id <- splitAsList(tx_id, sm)
    UATXHcount <- sgedges$UATXHcount
    if (!is.null(UATXHcount))
        sgedges$UATXHcount <- sum(splitAsList(sgedges0$UATXHcount, sm))
    if (is.null(ex_hits) && is.null(in_hits))
        return(sgedges)
    hits <- relist(character(0), PartitioningByEnd(NG=length(sm)))
    if (!is.null(ex_hits)) {
        if (!is(ex_hits, "CharacterList"))
            stop("'ex_hits' must be a CharacterList object")
        ex_idx <- which(ex_or_in == "ex")
        if (length(ex_idx) != length(ex_hits))
            stop("'ex_hits' is incompatible with 'sgedges0'")
        hits[ex_idx] <- ex_hits
    }
    if (!is.null(in_hits)) {
        if (!is(in_hits, "CharacterList"))
            stop("'in_hits' must be a CharacterList object")
        in_idx <- which(ex_or_in == "in")
        if (length(in_idx) != length(in_hits))
            stop("'in_hits' is incompatible with 'sgedges0'")
        hits[in_idx] <- in_hits
    }
    ## TODO: This is quite inefficient. Improve it.
    for (i in which(!is_not_dup))
        hits[[sm[i]]] <- unique(hits[[sm[i]]], hits[[i]])
    sgedges$hits <- hits[is_not_dup]
    sgedges$nhits <- elementLengths(sgedges$hits)
    sgedges
}

### Returns a DataFrame with first 2 cols being "from" and "to".
.extract_sgedges_exon_hits <- function(sg)
{
    if (length(sg) != 1L)
        stop("'sg' must be a SplicingGraphs object of length 1")
    ex_by_tx <- unlist(sg)
    exons <- ex_by_tx@unlistData
    common_strand <- commonStrand.GRanges(exons, what="exons in the gene")
    if (common_strand == "+") {
        from_to_colnames <- c("start_SSid", "end_SSid")
    } else {
        from_to_colnames <- c("end_SSid", "start_SSid")
    }
    ex_mcols <- mcols(exons)
    ex_colnames <- colnames(ex_mcols)
    hits_idx <- grep("hits$", ex_colnames)
    hits_colnames <- ex_colnames[hits_idx]
    hits_colnames <- c(from_to_colnames, hits_colnames)
    exon_hits <- ex_mcols[ , hits_colnames, drop=FALSE]
    colnames(exon_hits)[1:2] <- c("from", "to")
    exon_hits
}

### FIXME: Should return a DataFrame with first 2 cols being "from" and "to".
.extract_sgedges_intron_hits <- function(sg)
{
    if (length(sg) != 1L)
        stop("'sg' must be a SplicingGraphs object of length 1")
    in_by_tx <- intronsByTranscript(sg)
    introns <- in_by_tx@unlistData
    common_strand <- commonStrand.GRanges(introns, what="introns in the gene")
    ## FIXME: No "from" or "to" cols for now. Add them. This might require
    ## substantial changes upstream (i.e. SplicingGraphs() constructor) w.r.t
    ## what we store in the in_by_tx slot of the 'sg' object.
    #if (common_strand == "+") {
    #    from_to_colnames <- c("start_SSid", "end_SSid")
    #} else {
    #    from_to_colnames <- c("end_SSid", "start_SSid")
    #}
    in_mcols <- mcols(introns)
    in_colnames <- colnames(in_mcols)
    hits_idx <- grep("hits$", in_colnames)
    hits_colnames <- in_colnames[hits_idx]
    #hits_colnames <- c(from_to_colnames, hits_colnames)
    intron_hits <- in_mcols[ , hits_colnames, drop=FALSE]
    #colnames(intron_hits)[1:2] <- c("from", "to")
    intron_hits
}

setGeneric("sgedges", signature="x",
    function(x, UATXHcount=NULL, keep.dup.edges=FALSE)
        standardGeneric("sgedges")
)

setMethod("sgedges", "SplicingGraphs",
    function(x, UATXHcount=NULL, keep.dup.edges=FALSE)
    {
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        txpaths <- txpaths(x)
        if (is.null(UATXHcount))
            UATXHcount <- UATXHcount(x)
        if (keep.dup.edges)
            return(sgedges(txpaths, UATXHcount=UATXHcount,
                                    keep.dup.edges=keep.dup.edges))
        sgedges0 <- sgedges(txpaths, UATXHcount=UATXHcount,
                                     keep.dup.edges=TRUE)
        exon_hits <- .extract_sgedges_exon_hits(x)
        intron_hits <- .extract_sgedges_intron_hits(x)
        ## FIXME: Once .extract_sgedges_intron_hits() is fixed, merge the
        ## 'exon_hits' and 'intron_hits' DataFrame's with 'sgedges0' based
        ## on their "from" and "to" cols. Then remove the 'ex_hits' and
        ## 'in_hits' args from .make_sgedges_from_sgedges0().
        ex_hits <- exon_hits$hits
        in_hits <- intron_hits$hits
        .make_sgedges_from_sgedges0(sgedges0, ex_hits=ex_hits, in_hits=in_hits)
    }
)

setMethod("sgedges", "IntegerList",
    function(x, UATXHcount=NULL, keep.dup.edges=FALSE)
    {
        sgedges0 <- .make_sgedges0_from_txpaths(x, UATXHcount=UATXHcount)
        sgedges(sgedges0, keep.dup.edges=keep.dup.edges)
    }
)

setMethod("sgedges", "data.frame",
    function(x, UATXHcount=NULL, keep.dup.edges=FALSE)
    {
        if (!is.null(UATXHcount))
            stop("the 'UATXHcount' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        if (keep.dup.edges)
            return(x)  # no-op
        .make_sgedges_from_sgedges0(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgnodes() accessor
###

setGeneric("sgnodes", signature="x",
    function(x) standardGeneric("sgnodes")
)

setMethod("sgnodes", "ANY",
    function(x)
    {
        txpaths <- txpaths(x)
        sgnodes(txpaths)
    }
)

setMethod("sgnodes", "IntegerList",
    function(x) .get_sgnodes_from_txpaths(x)
)

setMethod("sgnodes", "data.frame",
    function(x) .get_sgnodes_from_sgedges(x)
)

setMethod("sgnodes", "DataFrame",
    function(x) .get_sgnodes_from_sgedges(x)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### outdeg() and indeg() extractors
###

setGeneric("outdeg", signature="x",
    function(x) standardGeneric("outdeg")
)

setMethod("outdeg", "ANY",
    function(x)
    {
        sgedges <- sgedges(x)
        outdeg(sgedges)
    }
)

setMethod("outdeg", "DataFrame",
    function(x)
    {
        sgnodes <- sgnodes(x)
        ans <- countMatches(sgnodes, x[ , "from"])
        names(ans) <- sgnodes
        ans
    }
)

setGeneric("indeg", signature="x",
    function(x) standardGeneric("indeg")
)

setMethod("indeg", "ANY",
    function(x)
    {
        sgedges <- sgedges(x)
        indeg(sgedges)
    }
)

setMethod("indeg", "DataFrame",
    function(x)
    {
        sgnodes <- sgnodes(x)
        ans <- countMatches(sgnodes, x[ , "to"])
        names(ans) <- sgnodes
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### uninformativeSSids() extractor
###
### Uninformative splicing sites are nodes in the splicing graph for which 
### outdeg and indeg are 1. A straightforward implementation of
### uninformativeSSids() would be:
###
###   uninformativeSSids <- function(x)
###   {
###       is_uninfo <- outdeg(sg) == 1L & indeg(sg) == 1L
###       names(is_uninfo)[is_uninfo]
###   }
###
### but the implementation below is about 2x faster.
###

setGeneric("uninformativeSSids", signature="x",
    function(x) standardGeneric("uninformativeSSids")
)

setMethod("uninformativeSSids", "ANY",
    function(x)
    {
        sgedges <- sgedges(x)
        uninformativeSSids(sgedges)
    }
)

setMethod("uninformativeSSids", "DataFrame",
    function(x)
    {
        from <- x[ , "from"]
        to <- x[ , "to"]
        from1_SSids <- setdiff(from, from[duplicated(from)])
        to1_SSids <- setdiff(to, to[duplicated(to)])
        intersect(from1_SSids, to1_SSids)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgedges2() extractor
###
### Same as sgedges() except that uninformative nodes (i.e. SSids) are removed.
###

### 'sgedges' must be a DataFrame as returned by:
###     sgedges( , keep.dup.edges=FALSE)
.remove_uninformative_SSids <- function(sgedges)
{
    ex_or_in <- sgedges[ , "ex_or_in"]
    ex_or_in_levels <- levels(ex_or_in)
    if (!identical(ex_or_in_levels, EX_OR_IN_LEVELS))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, the \"ex_or_in\" column has invalid ",
             "levels. Could it be that it was obtained by a previous call ",
             "to sgedges2()?")
    levels(ex_or_in) <- EX_OR_IN_LEVELS2
    uninformative_SSids <- uninformativeSSids(sgedges)
    if (length(uninformative_SSids) == 0L)
        return(sgedges)
    from <- sgedges[ , "from"]
    to <- sgedges[ , "to"]
    tx_id <- sgedges[ , "tx_id"]
    idx1 <- match(uninformative_SSids, from)
    idx2 <- match(uninformative_SSids, to)
    ## 2 sanity checks.
    if (!identical(unname(tx_id[idx1]), unname(tx_id[idx2])))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, the 2 rows containing a given ",
             "uninformative splicing site id must contain the same tx_id.",
             "Could it be that the \"tx_id\" column was manually altered ",
             "before the data.frame (or DataFrame) was passed to ",
             "sgedges2()?")
    if (!all(idx1 == idx2 + 1L))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, each uninformative splicing site ",
             "id must appear in 2 consecutive rows (first in the \"to\" ",
             "column, then in the \"from\" column. Could it be that the ",
             "rows were subsetted before the data.frame (or DataFrame) ",
             "was passed to sgedges2()?")
    from <- from[-idx1]
    to <- to[-idx2]
    ex_or_in[idx1] <- EX_OR_IN_LEVELS2[4L]
    ex_or_in <- ex_or_in[-idx2]
    tx_id <- tx_id[-idx1]
    DataFrame(from=from, to=to, ex_or_in=ex_or_in, tx_id=tx_id)
}

sgedges2 <- function(x)
{
    if (!is(x, "DataFrame"))
        x <- sgedges(x)
    .remove_uninformative_SSids(x)
}

