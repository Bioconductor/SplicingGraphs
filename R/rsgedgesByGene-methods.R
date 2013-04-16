### =========================================================================
### "rsgedgesByGene" (and related) methods
### -------------------------------------------------------------------------
###
### Analog to "sgedgesByGene" (and related) methods, but operate on the
### *reduced* splicing graph.
###

### Return the uninformative fully qualified nodes.
.uninformative_fqnodes <- function(sg)
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
    sgedge_id <- paste0(gene_id, ":", from, ",", to)
    keep_idx <- which(!duplicated(sgedge_id))
    fqfrom <- paste0(gene_id, ":", from)[keep_idx]  # fully qualified ids
    fqto <- paste0(gene_id, ":", to)[keep_idx]  # fully qualified ids
    uninformative_sgnodes(fqfrom, fqto)
}

.pmerge <- function(x, y)
{
    x_partitioning <- PartitioningByEnd(x)
    y_partitioning <- PartitioningByEnd(y)
    stopifnot(identical(x_partitioning, y_partitioning))
    starts <- start(x_partitioning)
    ends <- end(x_partitioning)
    unlisted_x <- unlist(x, use.names=FALSE)
    unlisted_y <- unlist(y, use.names=FALSE)
    stopifnot(identical(unlisted_x[-starts], unlisted_y[-ends]))
    y <- as(unlisted_y[ends], "List")
    xy <- c(x, y)
    f <- rep.int(seq_along(x), 2L)
    ans <- unlistAndSplit(xy, f)
    names(ans) <- names(x_partitioning)
    ans
}

### 'gene_id', 'from', and 'to', must have the same length N (nb of unique
### edges in the SplicingGraphs object befor reduction).
### Returns a character vector of length N containing the global rsgedge id
### (global reduced splicing graph edge id) corresponding to each input edge.
.build_sgedge2rsgedge_map <- function(gene_id, from, to, ui_fqnodes)
{
    ans <- paste0(gene_id, ":", from, ",", to)

    fqfrom <- paste0(gene_id, ":", from)  # fully qualified ids
    fqto <- paste0(gene_id, ":", to)  # fully qualified ids
    idx1 <- match(ui_fqnodes, fqfrom)
    idx2 <- match(ui_fqnodes, fqto)
    keep_idx <- which(!(is.na(idx1) | is.na(idx2)))
    idx1 <- idx1[keep_idx]
    idx2 <- idx2[keep_idx]
    stopifnot(all(idx1 == idx2 + 1L))  # sanity check
    group1 <- cumsum(c(TRUE, diff(idx1) != 1L))
    split_idx2 <- unname(splitAsList(idx2, group1))
    split_idx1 <- unname(splitAsList(idx1, group1))
    idx <- .pmerge(split_idx2, split_idx1)
    alter_idx <- idx@unlistData

    from_list <- relist(from[alter_idx], idx)
    to_list <- relist(to[alter_idx], idx)
    nodes_list <- .pmerge(from_list, to_list)

    rsgedge_id <- sapply(nodes_list, paste0, collapse=",")
    rsgedge_id <- rep.int(rsgedge_id, elementLengths(idx))
    rsgedge_id <- paste(gene_id[alter_idx], rsgedge_id, sep=":")

    ans[alter_idx] <- rsgedge_id
    ans
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

    ans_mcols <- DataFrame(rsgedge_id=levels(f),
                           from=ans_from,
                           to=ans_to,
                           ex_or_in=ans_ex_or_in,
                           tx_id=ans_tx_id)

    ## Reduce hits metadata cols.
    hits_mcol_idx <- grep("hits$", colnames(edges_mcols))
    if (length(hits_mcol_idx) != 0L) {
        ## FIXME: endoapply() on a DataFrame object is broken when applying
        ## a function 'FUN' that modifies the nb of rows. Furthermore, the
        ## returned object passes validitation despite being broken! Fix it
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
        ui_fqnodes <- .uninformative_fqnodes(x)
        sgedge2rsgedge_map <- .build_sgedge2rsgedge_map(gene_id, from, to,
                                                        ui_fqnodes)
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
             "to rsgedges()?")
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
             "rsgedges()?")
    if (!all(idx1 == idx2 + 1L))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, each uninformative splicing site ",
             "id must appear in 2 consecutive rows (first in the \"to\" ",
             "column, then in the \"from\" column. Could it be that the ",
             "rows were subsetted before the data.frame (or DataFrame) ",
             "was passed to rsgedges()?")
    from <- from[-idx1]
    to <- to[-idx2]
    ex_or_in[idx1] <- EX_OR_IN_LEVELS2[4L]
    ex_or_in <- ex_or_in[-idx2]
    tx_id <- tx_id[-idx1]
    DataFrame(from=from, to=to, ex_or_in=ex_or_in, tx_id=tx_id)
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

