### =========================================================================
### Summarizing the reads assigned to a SplicingGraphs object
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reportReads()
###

.reportReads_by_sgedge <- function(x)
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

.reportReads_by_rsgedge <- function(x)
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

.reportReads_by_tx <- function(x)
{
    ex_by_tx <- unlist(x)
    tx_id <- mcols(ex_by_tx)[ , "tx_id"]
    gene_id <- names(ex_by_tx)
    edges_by_tx <- sgedgesByTranscript(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_tx, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    hits_cols <- GenomicFeatures:::.collapse_df(edge_data[hits_mcol_idx],
                                                edges_by_tx)
    cbind(DataFrame(tx_id=tx_id, gene_id=gene_id), hits_cols)
}

.reportReads_by_gene <- function(x)
{
    edges_by_gene <- sgedgesByGene(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edges_by_gene, use.names=FALSE))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    hits_cols <- GenomicFeatures:::.collapse_df(edge_data[hits_mcol_idx],
                                                edges_by_gene)
    gene_id <- names(edges_by_gene)
    tx_id <- unique(IRanges:::regroupBySupergroup(edge_data[ , "tx_id"],
                                                  edges_by_gene))
    cbind(DataFrame(gene_id=gene_id, tx_id=tx_id), hits_cols)
}

setGeneric("reportReads", signature="x",
    function(x, by=c("sgedge", "rsgedge", "tx", "gene"))
        standardGeneric("reportReads")
)

### Return a DataFrame with 1 row per splicing graph edge (or reduced
### splicing graph edge), and 1 column per sample.
setMethod("reportReads", "SplicingGraphs",
    function(x, by=c("sgedge", "rsgedge", "tx", "gene"))
    {
        by <- match.arg(by)
        switch(by,
            sgedge=.reportReads_by_sgedge(x),
            rsgedge=.reportReads_by_rsgedge(x),
            tx=.reportReads_by_tx(x),
            gene=.reportReads_by_gene(x))
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
        reported_reads <- reportReads(x, by=by)
        hits_col_idx <- grep("\\.hits$", colnames(reported_reads))
        if (length(hits_col_idx) == 0L)
            return(reported_reads)
        read_counts <- endoapply(reported_reads[hits_col_idx], elementNROWS)
        colnames(read_counts) <- sub("\\.hits$", "", colnames(read_counts))
        cbind(reported_reads[-hits_col_idx], read_counts)
    }
)

