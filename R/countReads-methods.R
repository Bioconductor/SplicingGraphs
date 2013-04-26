### =========================================================================
### Summarizing the reads assigned to the edges of a SplicingGraphs object
### -------------------------------------------------------------------------


.countReads_by_sgedge <- function(x)
{
    edges_by_gene <- sgedgesByGene(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_gene, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    ans <- endoapply(edge_data[hits_mcol_idx], elementLengths)
    colnames(ans) <- sub("\\.hits$", "", colnames(ans))
    left_mcolnames <- c("sgedge_id", "ex_or_in")
    left_cols <- edge_data[left_mcolnames]
    cbind(left_cols, ans)
}

.countReads_by_rsgedge <- function(x)
{
    edges_by_gene <- rsgedgesByGene(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_gene, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    ans <- endoapply(edge_data[hits_mcol_idx], elementLengths)
    colnames(ans) <- sub("\\.hits$", "", colnames(ans))
    left_mcolnames <- c("rsgedge_id", "ex_or_in")
    left_cols <- edge_data[left_mcolnames]
    cbind(left_cols, ans)
}

.countReads_by_tx <- function(x)
{
    ex_by_tx <- unlist(x)
    tx_id <- mcols(ex_by_tx)[ , "tx_id"]
    gene_id <- names(ex_by_tx)
    edges_by_tx <- sgedgesByTranscript(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_tx, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    edge_data_breakpoints <- end(PartitioningByEnd(edges_by_tx))

    ## FIXME: endoapply() on a DataFrame object is broken when applying
    ## a function 'FUN' that modifies the nb of rows. Furthermore, the
    ## returned object passes validation despite being broken! Fix it
    ## in IRanges.
    ans <- endoapply(edge_data[hits_mcol_idx],
                     function(hits)
                       elementLengths(unique(regroup(hits,
                                                     edge_data_breakpoints))))

    ## Fix the broken DataFrame returned by endoapply().
    ans@nrows <- length(tx_id)
    ans@rownames <- NULL

    colnames(ans) <- sub("\\.hits$", "", colnames(ans))
    cbind(DataFrame(tx_id=tx_id, gene_id=gene_id), ans)
}

setGeneric("countReads", signature="x",
    function(x, by=c("sgedge", "rsgedge", "tx")) standardGeneric("countReads")
)

### Return a DataFrame with 1 row per splicing graph edge (or reduced
### splicing graph edge), and 1 column per sample.
setMethod("countReads", "SplicingGraphs",
    function(x, by=c("sgedge", "rsgedge", "tx"))
    {
        by <- match.arg(by)
        switch(by,
            sgedge=.countReads_by_sgedge(x),
            rsgedge=.countReads_by_rsgedge(x),
            tx=.countReads_by_tx(x))
    }
)

