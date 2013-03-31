### =========================================================================
### "sgedgesByTranscript" (and related) methods
### -------------------------------------------------------------------------


### Edge metadata columns that are considered to be exon attributes (note
### that we include the "start_SSid" and "end_SSid" cols). Those columns are
### the 5 first inner metadata columns of the GRangesList object containing
### the exons grouped by transcripts returned by unlist() when called on a
### SplicingGraphs object.
EXON_MCOLS <- c("exon_id", "exon_name", "exon_rank", "start_SSid", "end_SSid")

### All edge metadata columns.
ALL_EDGE_MCOLS <- c("sgedge_id", "from", "to", "ex_or_in", "tx_id", EXON_MCOLS)

### Subset of 'ALL_EDGE_MCOLS' made of those columns that are considered
### invariant i.e. the values in them associated with the same sgedge_id
### (global edge id) should be the same. Note that we also include the
### "sgedge_id" col itself.
INVARIANT_EDGE_MCOLS <- c("sgedge_id", "from", "to", "ex_or_in",
                          "start_SSid", "end_SSid")

EX_OR_IN_LEVELS2 <- c("ex", "in", "", "mixed")
EX_OR_IN_LEVELS <- EX_OR_IN_LEVELS2[-4L]

.check_exon_mcolnames <- function(colnames)
{
    stopifnot(identical(head(colnames, n=length(EXON_MCOLS)),
                        EXON_MCOLS))
}

.check_all_edge_mcolnames <- function(colnames)
{
    stopifnot(identical(head(colnames, n=length(ALL_EDGE_MCOLS)),
                        ALL_EDGE_MCOLS))
}

.get_index_of_mcols_to_remove <- function(colnames,
                                          with.exon.mcols, with.hits.mcols)
{
    ans <- integer(0)
    if (!with.exon.mcols) {
        idx <- match(EXON_MCOLS, colnames)
        ans <- c(ans, idx)
    }
    if (!with.hits.mcols) {
        idx <- grep("hits$", colnames)
        ans <- c(ans, idx)
    }
    ans
}

.get_index_of_invariant_edge_mcols <- function(colnames)
{
    idx <- match(INVARIANT_EDGE_MCOLS, colnames)
    idx[!is.na(idx)]
}


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
        .check_exon_mcolnames(colnames(ex_unlistData_mcols))

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
        .check_all_edge_mcolnames(colnames(ans_unlistData_mcols))

        ## Drop unwanted columns.
        mcol_idx <- .get_index_of_mcols_to_remove(
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

### .unlistAndSplit() example:
###
###   > x <- SimpleList(A=4:5, B=letters[1:4], C=NULL, D=1:2, E=-2:0, F=TRUE)
###   > f <- c("y", "x", "x", "x", "y", "x")
###   > .unlistAndSplit(x, f)
###   CharacterList of length 2
###   [["x"]] a b c d 1 2 TRUE
###   [["y"]] 4 5 -2 -1 0
###   > .unlistAndSplit(x, f)[[1]]
###        B      B      B      B      D      D      F 
###      "a"    "b"    "c"    "d"    "1"    "2" "TRUE" 
###
### Should work on any vector-like object and act as an endomorphism on a
### CompressedList object. On an atomic vector (on which 'unlist()' is a
### no-op), should be equivalent to 'splitAsList(x, f)'.
### TODO: Maybe move this to IRanges and expose to the user.
.unlistAndSplit <- function(x, f, drop=FALSE)
{
    if (length(f) != length(x))
        stop("'x' and 'f' must have the same length")
    x2 <- unlist(x)
    f2 <- rep.int(f, elementLengths(x))
    splitAsList(x2, f2, drop=drop)
}

setGeneric("sgedgesByGene", signature="x",
    function(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)
        standardGeneric("sgedgesByGene")
)

setMethod("sgedgesByGene", "SplicingGraphs",
    function(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)
    {
        edges_by_tx <- sgedgesByTranscript(x, with.exon.mcols=with.exon.mcols,
                                              with.hits.mcols=with.hits.mcols)
        edges0 <- unlist(edges_by_tx)
        edges0_mcols <- mcols(edges0)
        edges0_mcolnames <- colnames(edges0_mcols)

        sgedge_id <- edges0_mcols[ , "sgedge_id"]
        sm <- match(sgedge_id, sgedge_id)

        ## Sanity checks.
        stopifnot(all(edges0 == edges0[sm]))
        invariant_mcol_idx <- .get_index_of_invariant_edge_mcols(
                                  edges0_mcolnames)
        stopifnot(identical(
                    edges0_mcols[ , invariant_mcol_idx, drop=FALSE],
                    edges0_mcols[sm , invariant_mcol_idx, drop=FALSE]))

        ## Compute 'ans_partitioning'.
        keep_idx <- which(sm == seq_along(sm))
        ans_unlistData <- edges0[keep_idx]
        ans_grouping <- Rle(names(ans_unlistData))
        ans_eltlens <- runLength(ans_grouping)
        ans_partitioning <- PartitioningByEnd(cumsum(ans_eltlens),
                                              names=runValue(ans_grouping))

        ## Compute 'ans_unlistData'.
        names(ans_unlistData) <- NULL
        ans_unlistData_mcols <- mcols(ans_unlistData)

        variant_mcol_idx <- seq_along(edges0_mcolnames)[-invariant_mcol_idx]
        f <- factor(sgedge_id, levels=sgedge_id[keep_idx])
        for (i in variant_mcol_idx) {
            old_col <- edges0_mcols[ , i]
            new_col <- unname(unique(.unlistAndSplit(old_col, f)))
            ans_unlistData_mcols[ , i] <- new_col
        }
        mcols(ans_unlistData) <- ans_unlistData_mcols

        ## Relist 'ans_unlistData' and return.
        ans <- relist(ans_unlistData, ans_partitioning)
        ans
    }
)

