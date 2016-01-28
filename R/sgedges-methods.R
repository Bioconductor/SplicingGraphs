### =========================================================================
### "sgedges" (and related) methods
### -------------------------------------------------------------------------


.get_sgnodes_from_sgedges <- function(sgedges)
{
    from <- sgedges[ , "from"]
    to <- sgedges[ , "to"]
    SSids <- as.integer(setdiff(c(from, to), c("R", "L")))
    c("R", sort(SSids), "L")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgedges() extractor
###
### Returns the splicing graph in a DataFrame with 1 row per edge.
###

### 'txpath' must be an IntegerList containing all the splicing paths (1 per
### transcript) for a given gene. Should have been obtained thru the txpath()
### accessor. Returns a 4-col (or 5-col if 'txweight' is supplied) data.frame
### representing the splicing graph.
.make_sgedges0_from_txpath <- function(txpath, gene_id, txweight=NULL)
{
    if (!is.null(txweight)) {
        if (!is.numeric(txweight))
            stop("'txweight' must be a numeric vector or NULL")
        if (length(txweight) != length(txpath))
            stop("when not NULL, 'txweight' must have ",
                 "the same length as 'txpath'")
    }
    sgedges0s <- lapply(seq_along(txpath),
                        function(i) {
                            txpath_i <- txpath[[i]]
                            txpath_i_len <- length(txpath_i)
                            if (txpath_i_len %% 2L != 0L)
                                stop("some paths in 'txpath' contain ",
                                     "an odd number of splicing site ids")
                            from <- c("R", txpath_i)
                            to <- c(txpath_i, "L")
                            sgedge_id <- make_global_sgedge_id(gene_id,
                                                               from, to)
                            nexons <- txpath_i_len %/% 2L
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
                                       sgedge_id=sgedge_id,
                                       ex_or_in=ex_or_in,
                                       stringsAsFactors=FALSE)
                        })
    nedges_per_tx <- sapply(sgedges0s, nrow)
    sgedges0 <- do.call(rbind, sgedges0s)
    tx_id <- names(txpath)
    if (is.null(tx_id))
        tx_id <- seq_along(txpath)
    tx_id <- rep.int(tx_id, nedges_per_tx)
    sgedges0$tx_id <- tx_id
    if (!is.null(txweight))
        sgedges0$txweight <- rep.int(txweight, nedges_per_tx)
    sgedges0
}

### Collapse the duplicated edges in 'sgedges0' into a DataFrame.
### We use a DataFrame instead of a data.frame because we want to store
### the "tx_id" col in a CharacterList.
.make_sgedges_from_sgedges0 <- function(sgedges0, ex_hits=NULL, in_hits=NULL)
{
    from <- sgedges0[ , "from"]
    to <- sgedges0[ , "to"]
    sgedge_id <- sgedges0[ , "sgedge_id"]
    ex_or_in <- sgedges0[ , "ex_or_in"]
    tx_id <- sgedges0[ , "tx_id"]
    sm <- match(sgedge_id, sgedge_id)
    if (!all(ex_or_in == ex_or_in[sm]))
        stop("invalid splicing graph")
    is_not_dup <- sm == seq_along(sm)
    sgedges <- DataFrame(sgedges0[is_not_dup, , drop=FALSE])
    sgedges$tx_id <- unname(splitAsList(tx_id, sm))
    txweight <- sgedges$txweight
    if (!is.null(txweight))
        sgedges$txweight <- sum(splitAsList(sgedges0$txweight, sm))
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
    ## FIXME: Very inefficient. Improve it.
    for (i in which(!is_not_dup))
        hits[[sm[i]]] <- unique(hits[[sm[i]]], hits[[i]])
    hits <- hits[is_not_dup]
    sgedges$hits <- hits
    sgedges$nhits <- elementNROWS(hits)
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
    ex_mcolnames <- colnames(ex_mcols)
    hits_mcol_idx <- grep("hits$", ex_mcolnames)
    hits_mcolnames <- ex_mcolnames[hits_mcol_idx]
    hits_mcolnames <- c(from_to_colnames, hits_mcolnames)
    exon_hits <- ex_mcols[ , hits_mcolnames, drop=FALSE]
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
    hits_mcol_idx <- grep("hits$", in_colnames)
    hits_mcolnames <- in_colnames[hits_mcol_idx]
    #hits_mcolnames <- c(from_to_colnames, hits_mcolnames)
    intron_hits <- in_mcols[hits_mcolnames]
    #colnames(intron_hits)[1:2] <- c("from", "to")
    intron_hits
}

setGeneric("sgedges", signature="x",
    function(x, txweight=NULL, keep.dup.edges=FALSE)
        standardGeneric("sgedges")
)

setMethod("sgedges", "SplicingGraphs",
    function(x, txweight=NULL, keep.dup.edges=FALSE)
    {
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        txpath <- txpath(x)  # fails if length(x) != 1
        gene_id <- names(x)
        if (is.null(txweight))
            txweight <- txweight(x)
        sgedges0 <- .make_sgedges0_from_txpath(txpath, gene_id,
                                               txweight=txweight)
        if (keep.dup.edges)
            return(sgedges0)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgnodes() accessor
###

setGeneric("sgnodes", signature="x",
    function(x) standardGeneric("sgnodes")
)

setMethod("sgnodes", "SplicingGraphs",
    function(x)
    {
        txpath <- txpath(x)
        sgnodes(txpath)
    }
)

setMethod("sgnodes", "IntegerList",
    function(x) get_sgnodes_from_txpath(x)
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
        m <- match(x[ , "from"], sgnodes)
        ans <- tabulate(m, nbins=length(sgnodes))
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
        m <- match(x[ , "to"], sgnodes)
        ans <- tabulate(m, nbins=length(sgnodes))
        names(ans) <- sgnodes
        ans
    }
)

