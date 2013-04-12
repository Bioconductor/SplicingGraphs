### =========================================================================
### "sgedgesByGene" (and related) methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### intronsByTranscript()
###

setMethod("intronsByTranscript", "SplicingGraphs", function(x) x@in_by_tx)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgedgesByTranscript()
###

setGeneric("sgedgesByTranscript", signature="x",
    function(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)
        standardGeneric("sgedgesByTranscript")
)

setMethod("sgedgesByTranscript", "SplicingGraphs",
    function(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)
    {
        if (!isTRUEorFALSE(with.exon.mcols))
            stop("'with.exon.mcols' must be TRUE or FALSE")
        if (!isTRUEorFALSE(with.hits.mcols))
            stop("'with.hits.mcols' must be TRUE or FALSE")

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

        ## Add missing metadata cols to 'in_unlistData'.
        ex_unlistData <- ex_by_tx@unlistData
        ex_unlistData_len <- length(ex_unlistData)
        ex_unlistData_mcols <- mcols(ex_unlistData)
        check_exon_mcolnames(colnames(ex_unlistData_mcols))

        in_unlistData <- in_by_tx@unlistData
        in_unlistData_len <- length(in_unlistData)
        in_unlistData_mcols <- mcols(in_unlistData)

        in_missing_mcols <- DataFrame(exon_id=NA_integer_,
                                      exon_name=NA_character_,
                                      exon_rank=NA_integer_,
                                      start_SSid=NA_integer_,
                                      end_SSid=NA_integer_)
        idx <- rep.int(1L, in_unlistData_len)
        in_missing_mcols <- in_missing_mcols[idx, , drop=FALSE]
        in_unlistData_mcols <- cbind(in_missing_mcols, in_unlistData_mcols)

        ## Make "from" and "to" metadata cols.
        from <- as.character(ex_unlistData_mcols$start_SSid)
        to <- as.character(ex_unlistData_mcols$end_SSid)
        idx <- which(strand(ex_unlistData) == "-")
        tmp <- from[idx]
        from[idx] <- to[idx]
        to[idx] <- tmp
        ex_prepend_mcols <- DataFrame(from=from, to=to)

        from <- to <- rep.int(NA_character_, in_unlistData_len)
        in_prepend_mcols <- DataFrame(from=from, to=to)

        ## Make "ex_or_in" and "tx_id" metadata cols.
        ex_or_in <- rep.int(factor("ex", levels=EX_OR_IN_LEVELS),
                            ex_unlistData_len)
        ex_prepend_mcols$ex_or_in <- ex_or_in
        ex_or_in <- rep.int(factor("in", levels=EX_OR_IN_LEVELS),
                            in_unlistData_len)
        in_prepend_mcols$ex_or_in <- ex_or_in

        tx_id <- mcols(ex_by_tx)[ , "tx_id"]
        tx_id <- factor(tx_id, levels=unique(tx_id))
        ex_prepend_mcols$tx_id <- rep.int(tx_id, nex_by_tx)
        in_prepend_mcols$tx_id <- rep.int(tx_id, nin_by_tx)

        mcols(ex_unlistData) <- cbind(ex_prepend_mcols, ex_unlistData_mcols)
        mcols(in_unlistData) <- cbind(in_prepend_mcols, in_unlistData_mcols)

        ## Compute 'ans_unlistData'. We need to reorder 'c(ex_unlistData,
        ## in_unlistData)' to bring introns between their flanking exons.
        ans_unlistData <- c(ex_unlistData, in_unlistData)
        seq0 <- seq_along(ex_partitioning)
        roidx <- integer(ex_unlistData_len + in_unlistData_len)
        seq1 <- seq_len(ex_unlistData_len)
        idx1 <- seq1 * 2L - rep.int(seq0, nex_by_tx)
        roidx[idx1] <- seq1
        seq2 <- seq_len(in_unlistData_len)
        idx2 <- seq2 * 2L + rep.int(seq0, nin_by_tx) - 1L
        roidx[idx2] <- seq2 + ex_unlistData_len
        ans_unlistData <- ans_unlistData[roidx]

        ## Fill gaps in "from" and "to" metadata cols.
        ans_unlistData_mcols <- mcols(ans_unlistData)
        from <- ans_unlistData_mcols$from
        to <- ans_unlistData_mcols$to
        introns_idx <- which(is.na(from))  # same as 'which(is.na(to))'
        from[introns_idx] <- to[introns_idx - 1L]
        to[introns_idx] <- from[introns_idx + 1L]
        ans_unlistData_mcols$from <- from
        ans_unlistData_mcols$to <- to

        ## Sanity check: exons must be flanking introns.
        ans_unlistData_start <- start(ans_unlistData)
        ans_unlistData_end <- end(ans_unlistData)
        ans_unlistData_strand <- strand(ans_unlistData)

        plus_idx <- which(ans_unlistData_strand == "+")
        plus_introns_idx <- intersect(introns_idx, plus_idx)
        stopifnot(identical(ans_unlistData_start[plus_introns_idx] - 1L,
                            ans_unlistData_end[plus_introns_idx - 1L]))
        stopifnot(identical(ans_unlistData_end[plus_introns_idx] + 1L,
                            ans_unlistData_start[plus_introns_idx + 1L]))

        minus_idx <- which(ans_unlistData_strand == "-")
        minus_introns_idx <- intersect(introns_idx, minus_idx)
        stopifnot(identical(ans_unlistData_start[minus_introns_idx] - 1L,
                            ans_unlistData_end[minus_introns_idx + 1L]))
        stopifnot(identical(ans_unlistData_end[minus_introns_idx] + 1L,
                            ans_unlistData_start[minus_introns_idx - 1L]))

        ## Add "sgedge_id" metadata col.
        sgedge_id <- paste0(rep.int(gene_ids, width(ans_partitioning)), ":",
                            from, ",", to)
        ans_unlistData_mcols <- cbind(DataFrame(sgedge_id=sgedge_id),
                                      ans_unlistData_mcols)
        check_all_edge_mcolnames(colnames(ans_unlistData_mcols))

        ## Drop unwanted columns.
        mcol_idx <- get_index_of_group_of_mcols(
                        colnames(ans_unlistData_mcols),
                        with.exon.mcols, with.hits.mcols)
        if (length(mcol_idx) != 0L)
            ans_unlistData_mcols <- ans_unlistData_mcols[ , -mcol_idx,
                                                         drop=FALSE]
        mcols(ans_unlistData) <- ans_unlistData_mcols

        ## Relist 'ans_unlistData' and return.
        ans <- relist(ans_unlistData, ans_partitioning)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgedgesByGene()
###

setGeneric("sgedgesByGene", signature="x",
    function(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE,
                keep.dup.edges=FALSE)
        standardGeneric("sgedgesByGene")
)

setMethod("sgedgesByGene", "SplicingGraphs",
    function(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE,
                keep.dup.edges=FALSE)
    {
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        edges_by_tx <- sgedgesByTranscript(x, with.exon.mcols=with.exon.mcols,
                                              with.hits.mcols=with.hits.mcols)
        edges0 <- unlist(edges_by_tx)
        edges0_mcols <- mcols(edges0)
        edges0_mcolnames <- colnames(edges0_mcols)

        sgedge_id <- edges0_mcols[ , "sgedge_id"]
        sm <- match(sgedge_id, sgedge_id)

        ## Sanity checks.
        stopifnot(all(edges0 == edges0[sm]))
        invar_mcol_idx <- get_index_of_invariant_edge_mcols(edges0_mcolnames)
        stopifnot(identical(
                    edges0_mcols[ , invar_mcol_idx, drop=FALSE],
                    edges0_mcols[sm , invar_mcol_idx, drop=FALSE]))

        ## Compute 'ans_partitioning'.
        ans_unlistData <- edges0
        if (!keep.dup.edges) {
            keep_idx <- which(sm == seq_along(sm))
            ans_unlistData <- ans_unlistData[keep_idx]
        }
        ans_grouping <- Rle(names(ans_unlistData))
        ans_eltlens <- runLength(ans_grouping)
        ans_partitioning <- PartitioningByEnd(cumsum(ans_eltlens),
                                              names=runValue(ans_grouping))

        ## Compute 'ans_unlistData'.
        names(ans_unlistData) <- NULL
        if (!keep.dup.edges) {
            ans_unlistData_mcols <- mcols(ans_unlistData)
            var_mcol_idx <- seq_along(edges0_mcolnames)[- invar_mcol_idx]
            f <- factor(sgedge_id, levels=sgedge_id[keep_idx])
            for (i in var_mcol_idx) {
                old_col <- edges0_mcols[ , i]
                new_col <- unname(unique(unlistAndSplit(old_col, f)))
                ans_unlistData_mcols[ , i] <- new_col
            }
            mcols(ans_unlistData) <- ans_unlistData_mcols
        }

        ## Relist 'ans_unlistData' and return.
        ans <- relist(ans_unlistData, ans_partitioning)
        ans
    }
)

