### =========================================================================
### Summarizing the reads assigned to a SplicingGraphs object
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getReads()
###

.getReads_by_sgedge <- function(x)
{
    edges_by_gene <- sgedgesByGene(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_gene, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    hits_cols <- edge_data[hits_mcol_idx]
    left_mcolnames <- c("sgedge_id", "ex_or_in")
    left_cols <- edge_data[left_mcolnames]
    cbind(left_cols, hits_cols)
}

.getReads_by_rsgedge <- function(x)
{
    edges_by_gene <- rsgedgesByGene(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_gene, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    hits_cols <- edge_data[hits_mcol_idx]
    left_mcolnames <- c("rsgedge_id", "ex_or_in")
    left_cols <- edge_data[left_mcolnames]
    cbind(left_cols, hits_cols)
}

.getReads_by_tx <- function(x)
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
    hits_cols <- endoapply(edge_data[hits_mcol_idx],
                           function(hits)
                             unique(regroup(hits, edge_data_breakpoints)))

    ## Fix the broken DataFrame returned by endoapply().
    hits_cols@nrows <- length(tx_id)
    hits_cols@rownames <- NULL

    cbind(DataFrame(tx_id=tx_id, gene_id=gene_id), hits_cols)
}

.getReads_by_gene <- function(x)
{
    edges_by_gene <- sgedgesByGene(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_gene, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    edge_data_breakpoints <- end(PartitioningByEnd(edges_by_gene))

    ## FIXME: endoapply() on a DataFrame object is broken when applying
    ## a function 'FUN' that modifies the nb of rows. Furthermore, the
    ## returned object passes validation despite being broken! Fix it
    ## in IRanges.
    hits_cols <- endoapply(edge_data[hits_mcol_idx],
                     function(hits)
                       unique(regroup(hits, edge_data_breakpoints)))

    ## Fix the broken DataFrame returned by endoapply().
    hits_cols@nrows <- length(edges_by_gene)
    hits_cols@rownames <- NULL

    gene_id <- names(edges_by_gene)
    tx_id <- unique(regroup(edge_data[ , "tx_id"], edge_data_breakpoints))
    cbind(DataFrame(gene_id=gene_id, tx_id=tx_id), hits_cols)
}

setGeneric("getReads", signature="x",
    function(x, by=c("sgedge", "rsgedge", "tx", "gene"))
        standardGeneric("getReads")
)

### Return a DataFrame with 1 row per splicing graph edge (or reduced
### splicing graph edge), and 1 column per sample.
setMethod("getReads", "SplicingGraphs",
    function(x, by=c("sgedge", "rsgedge", "tx", "gene"))
    {
        by <- match.arg(by)
        switch(by,
            sgedge=.getReads_by_sgedge(x),
            rsgedge=.getReads_by_rsgedge(x),
            tx=.getReads_by_tx(x),
            gene=.getReads_by_gene(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countReads()
###

setGeneric("countReads", signature="x",
    function(x, by=c("sgedge", "rsgedge", "tx", "gene"))
        standardGeneric("countReads")
)

### Return a DataFrame with 1 row per splicing graph edge (or reduced
### splicing graph edge), and 1 column per sample.
setMethod("countReads", "SplicingGraphs",
    function(x, by=c("sgedge", "rsgedge", "tx", "gene"))
    {
        assigned_reads <- getReads(x, by=by)
        hits_col_idx <- grep("\\.hits$", colnames(assigned_reads))
        if (length(hits_col_idx) == 0L)
            return(assigned_reads)
        read_counts <- endoapply(assigned_reads[hits_col_idx], elementLengths)
        colnames(read_counts) <- sub("\\.hits$", "", colnames(read_counts))
        cbind(assigned_reads[-hits_col_idx], read_counts)
    }
)

