### =========================================================================
### "rsgedgesByGene" (and related) methods
### -------------------------------------------------------------------------
###
### Analog to "sgedgesByGene" (and related) methods, but operate on the
### *reduced* splicing graphs.
###

### Return the fully qualified uninformative nodes.
.get_fq_uninfnodes <- function(sg)
{
    txpath <- txpath(unlist(sg))
    skeleton2 <- PartitioningByEnd(end(PartitioningByEnd(txpath)) * 2L)
    txpath2 <- relist(rep(txpath@unlistData, each=2L), skeleton2)

    R_nodes <- rep.int(CharacterList("R"), length(txpath2))
    L_nodes <- rep.int(CharacterList("L"), length(txpath2))
    tmp <- c(R_nodes, as(txpath2, "CharacterList"), L_nodes)
    f <- rep.int(seq_along(txpath2), 3L)
    txpath3 <- unlistAndSplit(tmp, f)
    names(txpath3) <- names(txpath)

    gene_id <- rep.int(names(txpath3), elementLengths(txpath3))
    tmp <- txpath3@unlistData
    gene_id <- gene_id[c(TRUE, FALSE)]
    from <- tmp[c(TRUE, FALSE)]
    to <- tmp[c(FALSE, TRUE)]
    sgedge_id <- make_global_sgedge_id(gene_id, from, to)
    keep_idx <- which(!duplicated(sgedge_id))
    fq_from <- paste(gene_id, from, sep=":")[keep_idx]  # fully qualified ids
    fq_to <- paste(gene_id, to, sep=":")[keep_idx]  # fully qualified ids
    uninformative_sgnodes(fq_from, fq_to)
}

### 'f': factor. The "reverse factor" of 'f' is the named list of integer
### vectors that maps each level of 'f' to the positions in 'f' where that
### level is used. It can quickly be computed with:
###
###     revfactor <- split(seq_along(f), f)
###
### 'f' can be rebuilt from 'revfactor' with:
###
###     f2 <- .make_factor_from_revfactor(revfactor, length(f))
###     stopifnot(identical(f2, f)
###
.make_factor_from_revfactor <- function(revfactor, f_len)
{
    f_levels <- names(revfactor)
    idx <- integer(f_len)
    idx[] <- NA_integer_
    idx[unlist(revfactor, use.names=FALSE)] <- rep.int(seq_along(revfactor),
                                                   elementLengths(revfactor))
    factor(f_levels[idx], levels=f_levels)
}

### 'from' and 'to' must have the same length N (nb of unique edges in the
### SplicingGraphs object before reduction).
.make_revfactor_from_uninfnodes <- function(uninfnodes, from, to)
{
    from_idx <- match(uninfnodes, from)
    to_idx <- match(uninfnodes, to)
    keep_idx <- which(!(is.na(from_idx) | is.na(to_idx)))
    from_idx <- from_idx[keep_idx]
    to_idx <- to_idx[keep_idx]
    stopifnot(all(from_idx == to_idx + 1L))  # sanity check
    if (length(from_idx) == 0L) {
        group1 <- integer(0)
    } else {
        group1 <- cumsum(c(TRUE, diff(from_idx) != 1L))
    }
    split_to_idx <- unname(splitAsList(to_idx, group1))
    split_from_idx <- unname(splitAsList(from_idx, group1))
    fancy_punion(split_to_idx, split_from_idx)
}

### 'revfactor' must be a "reverse factor" as returned by
### .make_revfactor_from_uninfnodes().
### Returns a character vector of length N containing the global rsgedge id
### (global reduced splicing graph edge id) corresponding to each input edge.
.build_sgedge2rsgedge_map_from_revfactor <- function(revfactor,
                                                     gene_id, from, to)
{
    ans <- make_global_sgedge_id(gene_id, from, to)

    unlisted_revfactor <- unlist(revfactor, use.names=FALSE)
    from_list <- relist(from[unlisted_revfactor], revfactor)
    to_list <- relist(to[unlisted_revfactor], revfactor)
    nodes_list <- fancy_punion(from_list, to_list)

    rsgedge_id <- sapply(nodes_list, paste0, collapse=",")
    rsgedge_id <- rep.int(rsgedge_id, elementLengths(revfactor))
    rsgedge_id <- paste(gene_id[unlisted_revfactor], rsgedge_id, sep=":")

    ans[unlisted_revfactor] <- rsgedge_id
    ans
}

.build_sgedge2rsgedge_map <- function(uninfnodes, gene_id, from, to)
{
    fq_from <- paste(gene_id, from, sep=":")  # fully qualified ids
    fq_to <- paste(gene_id, to, sep=":")  # fully qualified ids
    revfactor <- .make_revfactor_from_uninfnodes(uninfnodes, fq_from, fq_to)
    .build_sgedge2rsgedge_map_from_revfactor(revfactor, gene_id, from, to)
}

### 'edges' must be a GRanges object.
.reduce_edges <- function(edges, f)
{
    if (!is(edges, "GRanges"))
        stop("'edges' must be a GRanges object")
    if (!is.factor(f))
        stop("'f' must be a factor")
    if (length(edges) != length(f))
        stop("'edges' and 'f' must have the same length")
    sm <- match(f, f)
    keep_idx <- sm == seq_along(sm)
    edges_names <- names(edges)
    edges_mcols <- mcols(edges)
    names(edges) <- mcols(edges) <- NULL
    edges_by_redge <- split(edges, f)
    partitioning <- PartitioningByEnd(edges_by_redge)
    partitioning_start <- start(partitioning)
    partitioning_end <- end(partitioning)

    ## The ranges within each top-level element of 'edges_by_redge' are an
    ## alternance of consecutive exons and introns from the same transcript
    ## and thus are adjacent (i.e. non-overlapping with no gaps between them).
    ## Therefore using reduce() or range() should be equivalent.
    redges_by_redge <- reduce(edges_by_redge)
    ans <- unlist(redges_by_redge)
    stopifnot(identical(names(ans), levels(f)))  # sanity check
    names(ans) <- edges_names[keep_idx]

    ## Reduce "from" and "to" metadata cols.
    edges_from <- edges_mcols[ , "from"]
    ans_from <- splitAsList(edges_from, f)
    edges_to <- edges_mcols[ , "to"]
    ans_to <- splitAsList(edges_to, f)
    stopifnot(identical(ans_from@unlistData[-partitioning_start],
                        ans_to@unlistData[-partitioning_end]))
    ans_from <- ans_from@unlistData[partitioning_start]
    ans_to <- ans_to@unlistData[partitioning_end]

    ## Reduce "ex_or_in" metadata col.
    edges_ex_or_in <- edges_mcols[ , "ex_or_in"]
    ans_ex_or_in <- splitAsList(edges_ex_or_in, f)
    levels(ans_ex_or_in@unlistData) <- EX_OR_IN_LEVELS2
    mixed_idx <- which(elementLengths(ans_ex_or_in) != 1L)
    ans_ex_or_in[mixed_idx] <- factor("mixed", levels=EX_OR_IN_LEVELS2)
    ans_ex_or_in <- factor(EX_OR_IN_LEVELS2[ans_ex_or_in@unlistData],
                           levels=EX_OR_IN_LEVELS2)

    ## Reduce "tx_id" metadata col.
    edges_tx_id <- edges_mcols[ , "tx_id"]
    stopifnot(identical(edges_tx_id, edges_tx_id[sm]))  # sanity check
    ans_tx_id <- edges_tx_id[keep_idx]

    ans_mcols <- DataFrame(from=ans_from,
                           to=ans_to,
                           rsgedge_id=levels(f),
                           ex_or_in=ans_ex_or_in,
                           tx_id=ans_tx_id)

    ## Reduce hits metadata cols.
    hits_mcol_idx <- grep("hits$", colnames(edges_mcols))
    if (length(hits_mcol_idx) != 0L) {
        ## FIXME: endoapply() on a DataFrame object is broken when applying
        ## a function 'FUN' that modifies the nb of rows. Furthermore, the
        ## returned object passes validation despite being broken! Fix it
        ## in IRanges.
        hits_mcols <- endoapply(edges_mcols[hits_mcol_idx],
                                function(hits)
                                  unname(unique(unlistAndSplit(hits, f))))

        ## Fix the broken DataFrame returned by endoapply().
        hits_mcols@nrows <- nlevels(f)
        hits_mcols@rownames <- NULL

        ## Combine with 'ans_mcols'.
        ans_mcols <- cbind(ans_mcols, hits_mcols)
    }

    ## Set the metadata cols on 'ans' and return it.
    mcols(ans) <- ans_mcols
    ans
}


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
        uninformative_sgnodes(from, to)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rsgedgesByTranscript()
###

setGeneric("rsgedgesByTranscript", signature="x",
    function(x, with.hits.mcols=FALSE)
        standardGeneric("rsgedgesByTranscript")
)

setMethod("rsgedgesByTranscript", "SplicingGraphs",
    function(x, with.hits.mcols=FALSE)
    {
        stop("not ready yet!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rsgedgesByGene()
###

setGeneric("rsgedgesByGene", signature="x",
    function(x, with.hits.mcols=FALSE, keep.dup.edges=FALSE)
        standardGeneric("rsgedgesByGene")
)

setMethod("rsgedgesByGene", "SplicingGraphs",
    function(x, with.hits.mcols=FALSE, keep.dup.edges=FALSE)
    {
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        if (keep.dup.edges)
            stop("'keep.dup.edges=TRUE' is not supported yet, sorry")
        edges_by_gene <- sgedgesByGene(x, with.hits.mcols=with.hits.mcols)
        edges0 <- unlist(edges_by_gene)
        edges0_mcols <- mcols(edges0)
        gene_id <- names(edges0)
        from <- edges0_mcols[ , "from"]
        to <- edges0_mcols[ , "to"]
        uninfnodes <- .get_fq_uninfnodes(x)
        sgedge2rsgedge_map <- .build_sgedge2rsgedge_map(uninfnodes,
                                                        gene_id, from, to)
        f <- factor(sgedge2rsgedge_map, levels=unique(sgedge2rsgedge_map))
        ans_flesh <- .reduce_edges(edges0, f)
        ans_flesh_names <- Rle(names(ans_flesh))
        breakpoints <- cumsum(runLength(ans_flesh_names))
        ans_names <- runValue(ans_flesh_names)
        ans_skeleton <- PartitioningByEnd(breakpoints, names=ans_names)
        relist(unname(ans_flesh), ans_skeleton)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rsgedges() extractor
###
### Same as sgedges() except that uninformative nodes are removed.
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
             "to rsgedges()?")
    levels(ex_or_in) <- EX_OR_IN_LEVELS2
    uninformative_SSids <- uninformativeSSids(sgedges)
    if (length(uninformative_SSids) == 0L) {
        col_idx <- match("sgedge_id", colnames(sgedges))
        colnames(sgedges)[col_idx] <- "rsgedge_id"
        sgedges$ex_or_in <- ex_or_in
        return(sgedges)
    }

    from <- sgedges[ , "from"]
    to <- sgedges[ , "to"]
    from_idx <- match(uninformative_SSids, from)
    to_idx <- match(uninformative_SSids, to)
    if (!all(from_idx == to_idx + 1L))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, each uninformative splicing site ",
             "id must appear in 2 consecutive rows (first in the \"to\" ",
             "column, then in the \"from\" column. Could it be that the ",
             "rows were subsetted before the data.frame (or DataFrame) ",
             "was passed to rsgedges()?")

    ## Reduce "from" and "to" cols.
    ans_from <- from[-from_idx]
    ans_to <- to[-to_idx]

    ## Reduce "sgedge_id" col.
    sgedges_id <- sgedges[ , "sgedge_id"]
    tmp <- unlist(strsplit(sgedges_id, ":", fixed=TRUE), use.names=FALSE)
    gene_id <- tmp[c(TRUE, FALSE)]
    revfactor <- .make_revfactor_from_uninfnodes(uninformative_SSids, from, to)
    sgedge2rsgedge_map <- .build_sgedge2rsgedge_map_from_revfactor(revfactor,
                                                             gene_id, from, to)
    f <- factor(sgedge2rsgedge_map, levels=unique(sgedge2rsgedge_map))
    ans_rsgedge_id <- levels(f)

    ## Reduce "ex_or_in" col.
    ex_or_in[from_idx] <- EX_OR_IN_LEVELS2[4L]
    ans_ex_or_in <- ex_or_in[-to_idx]

    ## Reduce "tx_id" col.
    tx_id <- sgedges[ , "tx_id"]
    if (!identical(unname(tx_id[from_idx]), unname(tx_id[to_idx])))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, the 2 rows containing a given ",
             "uninformative splicing site id must contain the same tx_id.",
             "Could it be that the \"tx_id\" column was manually altered ",
             "before the data.frame (or DataFrame) was passed to ",
             "rsgedges()?")
    ans_tx_id <- tx_id[-from_idx]

    DataFrame(from=ans_from,
              to=ans_to,
              rsgedge_id=ans_rsgedge_id,
              ex_or_in=ans_ex_or_in,
              tx_id=ans_tx_id)
}

rsgedges <- function(x)
{
    if (!is(x, "DataFrame"))
        x <- sgedges(x)
    .remove_uninformative_SSids(x)
}

### Alias for backward compatibility.
sgedges2 <- rsgedges


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rsgraph() extractor
###
### Same as sgraph() except that uninformative nodes (i.e. SSids) are removed.
###

rsgraph <- function(x, tx_id.as.edge.label=FALSE, as.igraph=FALSE)
{
    sgraph(rsgedges(x),
           tx_id.as.edge.label=tx_id.as.edge.label, as.igraph=as.igraph)
}

### Alias for backward compatibility.
sgraph2 <- rsgraph

