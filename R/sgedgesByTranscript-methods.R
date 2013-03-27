### =========================================================================
### "sgedgesByTranscript" methods
### -------------------------------------------------------------------------


setGeneric("sgedgesByTranscript",
    function(x) standardGeneric("sgedgesByTranscript")
)

setMethod("sgedgesByTranscript", "SplicingGraphs",
    function(x)
    {
        ex_by_tx <- unlist(x)
        ex_partitioning <- PartitioningByEnd(ex_by_tx)
        gene_ids <- names(ex_partitioning)
        in_by_tx <- intronsByTranscript(x)
        in_partitioning <- PartitioningByEnd(in_by_tx)
        stopifnot(identical(gene_ids, names(in_partitioning)))

        ## Compute 'ans_partitioning'.
        nex_by_tx <- width(ex_partitioning)
        nin_by_tx <- width(in_partitioning)
        stopifnot(identical(nin_by_tx + 1L, nex_by_tx))
        ans_partitioning <- PartitioningByEnd(end(ex_partitioning) +
                                              end(in_partitioning),
                                              names=gene_ids)

        ## Compute 'ans_unlistData'.
        ex_unlistData <- ex_by_tx@unlistData
        ex_unlistData_len <- length(ex_unlistData)
        in_unlistData <- in_by_tx@unlistData
        in_unlistData_len <- length(in_unlistData)
        in_missing_mcols <- DataFrame(exon_id=NA_integer_,
                                      exon_name=NA_character_,
                                      exon_rank=NA_integer_,
                                      start_SSid=NA_integer_,
                                      end_SSid=NA_integer_)
        idx <- rep.int(1L, in_unlistData_len)
        in_missing_mcols <- in_missing_mcols[idx, , drop=FALSE]
        mcols(in_unlistData) <- cbind(in_missing_mcols, mcols(in_unlistData))
        ans_unlistData <- c(ex_unlistData, in_unlistData)

        ## Reorder 'ans_unlistData' to bring introns between the corresponding
        ## exons.
        seq0 <- seq_along(ex_partitioning)
        roidx <- integer(ex_unlistData_len + in_unlistData_len)
        seq1 <- seq_len(ex_unlistData_len)
        idx1 <- seq1 * 2L - rep.int(seq0, nex_by_tx)
        roidx[idx1] <- seq1
        seq2 <- seq_len(in_unlistData_len)
        idx2 <- seq2 * 2L + rep.int(seq0, nin_by_tx) - 1L
        roidx[idx2] <- seq2 + ex_unlistData_len
        ans_unlistData <- ans_unlistData[roidx]

        ## Relist 'ans_unlistData' and return.
        ans <- relist(ans_unlistData, ans_partitioning)
        ans
    }
)

