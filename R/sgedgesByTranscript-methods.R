### =========================================================================
### "sgedgesByTranscript" methods
### -------------------------------------------------------------------------


EX_OR_IN_LEVELS2 <- c("ex", "in", "", "mixed")
EX_OR_IN_LEVELS <- EX_OR_IN_LEVELS2[-4L]


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

        ## Add missing metadata cols to 'in_unlistData'.
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

        ## Add "from" and "to" metadata cols to 'ex_unlistData'
        ## and 'in_unlistData'.
        ex_unlistData_mcols <- mcols(ex_unlistData)
        from <- as.character(ex_unlistData_mcols$start_SSid)
        to <- as.character(ex_unlistData_mcols$end_SSid)
        idx <- which(strand(ex_unlistData) == "-")
        tmp <- from[idx]
        from[idx] <- to[idx]
        to[idx] <- tmp
        mcols(ex_unlistData)$from <- from
        mcols(ex_unlistData)$to <- to

        from <- to <- rep.int(NA_character_, in_unlistData_len)
        mcols(in_unlistData)$from <- from
        mcols(in_unlistData)$to <- to

        ## Add "ex_or_in" and "tx_id" metadata cols to 'ex_unlistData'
        ## and 'in_unlistData'.
        mcols(ex_unlistData)$ex_or_in <- factor("ex", levels=EX_OR_IN_LEVELS)
        mcols(in_unlistData)$ex_or_in <- factor("in", levels=EX_OR_IN_LEVELS)
        tx_id <- mcols(ex_by_tx)[ , "tx_id"]
        tx_id <- factor(tx_id, levels=unique(tx_id))
        mcols(ex_unlistData)$tx_id <- rep.int(tx_id, nex_by_tx)
        mcols(in_unlistData)$tx_id <- rep.int(tx_id, nin_by_tx)

        ## Compute 'ans_unlistData'. We need to reorder 'c(ex_unlistData,
        ## in_unlistData)' to bring introns between the corresponding exons.
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
        mcols(ans_unlistData)$from <- from
        mcols(ans_unlistData)$to <- to

        ## Sanity check: comparing intron starts/ends with exon starts/ends.
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

        ## Relist 'ans_unlistData' and return.
        ans <- relist(ans_unlistData, ans_partitioning)
        ans
    }
)

