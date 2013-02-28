### =========================================================================
### sgdf (and related) methods
### -------------------------------------------------------------------------


EX_OR_IN_LEVELS2 <- c("ex", "in", "", "mixed")
EX_OR_IN_LEVELS <- EX_OR_IN_LEVELS2[-4L]


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spath() accessor
###
### Gets all the splicing paths for the specified gene.
### Returns them in a named IntegerList with 1 top-level element per
### transcript in the specified gene. Each top-level element 'spath[[i]]'
### contains the splicing site ids for the i-th transcript.
###

setGeneric("spath", signature="x",
    function(x, gene_id=NA) standardGeneric("spath")
)

### Should return a CompressedIntegerList.
setMethod("spath", "SplicingGraphs",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        ans <- mcols(x@tx)[ , "spath"]
        if (is.null(x_names)) {
            if (!is.na(gene_id))
                stop("the 'gene_id' arg is not supported ",
                     "when 'x' is unnamed (in which case all its elements ",
                     "(i.e. transcripts) are considered to belong to the ",
                     "same gene)")
            return(ans) 
        }
        if (is.na(gene_id))
            stop("'gene_id' must be supplied when 'x' has names")
        ans <- ans[x_names == gene_id]
        if (length(ans) == 0L)
            stop("invalid 'gene_id'")
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UATXHcount() accessor
###

setGeneric("UATXHcount", signature="x",
    function(x, gene_id=NA) standardGeneric("UATXHcount")
)

### Should return an integer vector or a NULL.
setMethod("UATXHcount", "SplicingGraphs",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        ans <- mcols(x@tx)[["UATXHcount"]]
        if (is.null(x_names)) {
            if (!is.na(gene_id))
                stop("the 'gene_id' arg is not supported ",
                     "when 'x' is unnamed (in which case all its elements ",
                     "(i.e. transcripts) are considered to belong to the ",
                     "same gene)")
            return(ans) 
        }
        if (is.na(gene_id))
            stop("'gene_id' must be supplied when 'x' has names")
        if (is.null(ans))
            return(ans)
        ans <- ans[x_names == gene_id]
        if (length(ans) == 0L)
            stop("invalid 'gene_id'")
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .hits() accessor (not exported)
###

setGeneric(".hits", signature="x",
    function(x, gene_id=NA) standardGeneric(".hits")
)

### Should return a CompressedCharacterList or a NULL.
setMethod(".hits", "GRangesList",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        if (is.null(x_names)) {
            if (!is.na(gene_id))
                stop("the 'gene_id' arg is not supported ",
                     "when 'x' is unnamed (in which case all its elements ",
                     "(i.e. transcripts) are considered to belong to the ",
                     "same gene)")
            ans <- mcols(unlist(x, use.names=FALSE))[["hits"]]
            return(ans) 
        }
        if (is.na(gene_id))
            stop("'gene_id' must be supplied when 'x' has names")
        x <- x[x_names == gene_id]
        if (length(x) == 0L)
            stop("invalid 'gene_id'")
        ans <- mcols(unlist(x, use.names=FALSE))[["hits"]]
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgdf() extractor
###
### Returns the splicing graph in a DataFrame with 1 row per edge.
###

### 'spath' must be an IntegerList containing all the splicing paths for a
### given gene. Should have been obtained thru the spath() accessor.
### Returns a 4-col (or 5-col if 'UATXHcount' is supplied) data.frame
### representing the splicing graph.
.make_sgdf0_from_spath <- function(spath, UATXHcount=NULL)
{
    if (!is.null(UATXHcount)) {
        if (!is.integer(UATXHcount))
            stop("'UATXHcount' must be an integer vector or NULL")
        if (length(UATXHcount) != length(spath))
            stop("when not NULL, 'UATXHcount' must have ",
                 "the same length as 'spath'")
    }
    sgdf0s <- lapply(seq_along(spath),
                     function(i) {
                         SSids <- spath[[i]]
                         from <- c("R", SSids)
                         to <- c(SSids, "L")
                         nb_SSids <- length(SSids)
                         if (nb_SSids %% 2L != 0L)
                             stop("some splicing paths in 'spath' go thru an ",
                                  "odd number of splicing site ids")
                         nexons <- nb_SSids %/% 2L
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
    nedges_per_tx <- sapply(sgdf0s, nrow)
    sgdf0 <- do.call(rbind, sgdf0s)
    tx_id <- names(spath)
    if (is.null(tx_id))
        tx_id <- seq_along(spath)
    tx_id <- rep.int(factor(tx_id, levels=tx_id), nedges_per_tx)
    sgdf0$tx_id <- tx_id
    if (!is.null(UATXHcount))
        sgdf0$UATXHcount <- rep.int(UATXHcount, nedges_per_tx)
    sgdf0
}

### Collapse the duplicated edges in 'sgdf0' into a DataFrame.
### We use a DataFrame instead of a data.frame because we want to store
### the tx_id col in a CompressedFactorList (even though this container
### doesn't formally exist and a CompressedIntegerList is actually what's
### being used).
.make_sgdf_from_sgdf0 <- function(sgdf0, ex_hits=NULL, in_hits=NULL)
{
    from <- sgdf0[ , "from"]
    to <- sgdf0[ , "to"]
    ex_or_in <- sgdf0[ , "ex_or_in"]
    tx_id <- sgdf0[ , "tx_id"]
    edges <- paste(from, to, sep="~")
    sm <- match(edges, edges)
    if (!all(ex_or_in == ex_or_in[sm]))
        stop("invalid splicing graph")
    is_not_dup <- sm == seq_along(sm)
    sgdf <- DataFrame(sgdf0[is_not_dup, , drop=FALSE])
    sgdf$tx_id <- splitAsList(tx_id, sm)
    UATXHcount <- sgdf$UATXHcount
    if (!is.null(UATXHcount))
        sgdf$UATXHcount <- sum(splitAsList(sgdf0$UATXHcount, sm))
    if (is.null(ex_hits) && is.null(in_hits))
        return(sgdf)
    hits <- relist(character(0), PartitioningByEnd(NG=length(sm)))
    if (!is.null(ex_hits)) {
        if (!is(ex_hits, "CharacterList"))
            stop("'ex_hits' must be a CharacterList object")
        ex_idx <- which(ex_or_in == "ex")
        if (length(ex_idx) != length(ex_hits))
            stop("'ex_hits' is incompatible with 'sgdf0'")
        hits[ex_idx] <- ex_hits
    }
    if (!is.null(in_hits)) {
        if (!is(in_hits, "CharacterList"))
            stop("'in_hits' must be a CharacterList object")
        in_idx <- which(ex_or_in == "in")
        if (length(in_idx) != length(in_hits))
            stop("'in_hits' is incompatible with 'sgdf0'")
        hits[in_idx] <- in_hits
    }
    ## TODO: This is quite inefficient. Improve it.
    for (i in which(!is_not_dup))
        hits[[sm[i]]] <- unique(hits[[sm[i]]], hits[[i]])
    sgdf$hits <- hits[is_not_dup]
    sgdf$nhits <- elementLengths(sgdf$hits)
    sgdf
}

setGeneric("sgdf", signature="x",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
        standardGeneric("sgdf")
)

setMethod("sgdf", "ANY",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
    {
        spath <- spath(x, gene_id=gene_id)
        if (is.null(UATXHcount))
            UATXHcount <- UATXHcount(x, gene_id=gene_id)
        if (is.null(inbytx))
            return(sgdf(spath, UATXHcount=UATXHcount,
                               keep.dup.edges=keep.dup.edges))
        if (!is(inbytx, "GRangesList"))
            stop("'inbytx' must be NULL or a GRangesList object")
        if (!is(x, "SplicingGraphs"))
            stop("'x' must be a SplicingGraphs object ",
                 "when 'inbytx' is a GRangesList object")
        if (length(inbytx) != length(x))
            stop("'inbytx' must have the same length as 'x'")
        if (!identical(elementLengths(inbytx) + 1L, elementLengths(x)))
            stop("the shape of 'inbytx' is not compatible ",
                 "with the shape of 'x'")
        if (!identical(keep.dup.edges, FALSE))
            stop("'keep.dup.edges' must be FALSE when 'inbytx' is supplied")
        sgdf0 <- sgdf(spath, UATXHcount=UATXHcount, keep.dup.edges=TRUE)
        ex_or_in <- sgdf0[ , "ex_or_in"]
        ex_hits <- .hits(x@tx, gene_id=gene_id)
        if (is.null(ex_hits))
            stop("'x' must have a \"hits\" inner metadata column ",
                 "when 'inbytx' is a GRangesList object. May be ",
                 "you forgot to pass it thru assignSubfeatureHits()?")
        in_hits <- .hits(inbytx, gene_id=gene_id)
        if (is.null(in_hits))
            stop("'inbytx' has no \"hits\" inner metadata column. May be ",
                 "you forgot to pass it thru assignSubfeatureHits()?")
        .make_sgdf_from_sgdf0(sgdf0, ex_hits=ex_hits, in_hits=in_hits)
    }
)

setMethod("sgdf", "IntegerList",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an IntegerList")
        if (!is.null(inbytx))
            stop("the 'inbytx' arg is not supported ",
                 "when 'x' is an IntegerList")
        sgdf0 <- .make_sgdf0_from_spath(x, UATXHcount=UATXHcount)
        sgdf(sgdf0, keep.dup.edges=keep.dup.edges)
    }
)

setMethod("sgdf", "data.frame",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!is.null(UATXHcount))
            stop("the 'UATXHcount' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!is.null(inbytx))
            stop("the 'inbytx' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        if (keep.dup.edges)
            return(x)  # no-op
        .make_sgdf_from_sgdf0(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### uninformativeSSids() extractor
###

setGeneric("uninformativeSSids", signature="x",
    function(x, gene_id=NA) standardGeneric("uninformativeSSids")
)

setMethod("uninformativeSSids", "ANY",
    function(x, gene_id=NA)
    {
        x <- sgdf(x, gene_id=gene_id)
        uninformativeSSids(x)
    }
)

setMethod("uninformativeSSids", "DataFrame",
    function(x, gene_id=NA)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a DataFrame")
        from <- x[ , "from"]
        to <- x[ , "to"]
        from1_SSids <- setdiff(from, from[duplicated(from)])
        to1_SSids <- setdiff(to, to[duplicated(to)])
        intersect(from1_SSids, to1_SSids)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgdf2() extractor
###
### Same as sgdf() except that uninformative nodes (i.e. SSids) are removed.
###

### 'sgdf' must be a DataFrame as returned by:
###     sgdf( , keep.dup.edges=FALSE)
.remove_uninformative_SSids <- function(sgdf)
{
    ex_or_in <- sgdf[ , "ex_or_in"]
    ex_or_in_levels <- levels(ex_or_in)
    if (!identical(ex_or_in_levels, EX_OR_IN_LEVELS))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, the \"ex_or_in\" column has invalid ",
             "levels. Could it be that it was obtained by a previous call ",
             "to sgdf2()?")
    levels(ex_or_in) <- EX_OR_IN_LEVELS2
    uninformative_SSids <- uninformativeSSids(sgdf)
    if (length(uninformative_SSids) == 0L)
        return(sgdf)
    from <- sgdf[ , "from"]
    to <- sgdf[ , "to"]
    tx_id <- sgdf[ , "tx_id"]
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
             "sgdf2()?")
    if (!all(idx1 == idx2 + 1L))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, each uninformative splicing site ",
             "id must appear in 2 consecutive rows (first in the \"to\" ",
             "column, then in the \"from\" column. Could it be that the ",
             "rows were subsetted before the data.frame (or DataFrame) ",
             "was passed to sgdf2()?")
    from <- from[-idx1]
    to <- to[-idx2]
    ex_or_in[idx1] <- EX_OR_IN_LEVELS2[4L]
    ex_or_in <- ex_or_in[-idx2]
    tx_id <- tx_id[-idx1]
    DataFrame(from=from, to=to, ex_or_in=ex_or_in, tx_id=tx_id)
}

sgdf2 <- function(x, gene_id=NA)
{
    if (!is(x, "DataFrame"))
        x <- sgdf(x, gene_id=gene_id)
    .remove_uninformative_SSids(x)
}

