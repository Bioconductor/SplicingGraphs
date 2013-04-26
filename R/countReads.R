### =========================================================================
### Functions for assigning reads to the edges of a SplicingGraphs object and
### for summarizing them
### -------------------------------------------------------------------------


.check_reads_names <- function(reads)
{
    reads_names <- names(reads)
    if (is.null(reads_names))
        stop("'reads' must have names")
    errmsg <- c("'reads' has duplicated names. This probably means that",
                "either: (a) it contains paired-end reads that have not",
                "been paired, or (b) some of the reads are secondary",
                "alignments. See 'About the read names' section in the man",
                "page for assignReads() (accessible with '?assignReads')",
                "for more information.")
    if (anyDuplicated(reads_names))
        stop(paste0(errmsg, collapse="\n  "))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .assignSubfeatureHits()
###
### This is the workhorse behind assignReads().
###

### 'reads': a named GRangesList object containing gapped reads.
###     Alternatively it can be the GRanges object obtained by extracting
###     a single range per read e.g. the spanning ranges obtained with
###     'unlist(range(reads))'. In any case, it must have the length and
###     names of the original GRangesList object.
### 'subject': a GRangesList object containing the ranges of some subfeature
###     (e.g. exonic or intronic ranges) grouped by their parent feature
###     (e.g. transcripts).
### 'hits': a Hits object compatible with 'reads' and 'subject'.
.check_assignSubfeatureHits_args <- function(reads, subject, hits,
                                             ignore.strand, hits.colname)
{
    if (!(is(reads, "GRangesList") || is(reads, "GRanges")))
        stop("'reads' must be a GRangesList or GRanges object")
    .check_reads_names(reads)
    if (!is(subject, "GRangesList"))
        stop("'subject' must be a GRangesList object")
    if (!is(hits, "Hits"))
        stop("'hits' must be a Hits object")
    if (queryLength(hits) != length(reads)
     || subjectLength(hits) !=  length(subject))
        stop("'hits' is not compatible with 'reads' and 'subject'")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    if (!isSingleString(hits.colname))
        stop("'hits.colname' must be a single string")
}

### Returns 'subject' with 1 additional inner metadata col "hits" containing
### the hits assigned to each subrange.
.assignSubfeatureHits <- function(reads, subject, hits, ignore.strand=FALSE,
                                  hits.colname="hits")
{
    .check_assignSubfeatureHits_args(reads, subject, hits,
                                     ignore.strand, hits.colname)
    reads_names <- names(reads)
    subject_unlistData <- subject@unlistData
    subhits <- findOverlaps(reads, subject_unlistData,
                            ignore.strand=ignore.strand)
    subhits_q <- queryHits(subhits)
    subhits_s <- togroup(subject@partitioning, subjectHits(subhits))
    m <- IRanges:::matchIntegerPairs(subhits_q, subhits_s,
                                     queryHits(hits), subjectHits(hits))
    subhits <- subhits[!is.na(m)]
    hit_per_subfeature <- splitAsList(reads_names[queryHits(subhits)],
                                      subjectHits(subhits))
    mcols(subject_unlistData)[[hits.colname]] <- CharacterList(character(0))
    idx <- as.integer(names(hit_per_subfeature))
    names(hit_per_subfeature) <- NULL
    mcols(subject_unlistData)[[hits.colname]][idx] <- hit_per_subfeature
    subject@unlistData <- subject_unlistData
    subject
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### assignReads()
###

assignReads <- function(sg, reads, sample.name=NA)
{
    if (!is(sg, "SplicingGraphs"))
        stop("'sg' must be a SplicingGraphs object")
    if (is(reads, "GAlignments") || is(reads, "GAlignmentPairs")) {
        reads <- grglist(reads, order.as.in.query=TRUE)
    } else if (!is(reads, "GRangesList")) {
        stop("'reads' must be a GAlignments, GAlignmentPairs, ",
             "or GRangesList object")
    }
    if (!isSingleStringOrNA(sample.name))
        stop("'sample.name' must be a single string or NA")
    if (is.na(sample.name)) {
        hits.colname <- "hits"
    } else {
        hits.colname <- paste0(sample.name, ".hits")
    }

    ex_by_tx <- sg@genes@unlistData
    ov0 <- findOverlaps(reads, ex_by_tx, ignore.strand=TRUE)
    ovenc0 <- encodeOverlaps(reads, ex_by_tx, hits=ov0,
                             flip.query.if.wrong.strand=TRUE)
    ov0_is_compat <- isCompatibleWithSplicing(ovenc0)
    ov1 <- ov0[ov0_is_compat]

    query.breaks <- mcols(reads)$query.break
    ## 'reads2' is obtained by removing the gaps (i.e. Ns in the CIGAR) in
    ## 'reads'.
    if (is.null(query.breaks)) {
        ## Single-end reads. We produce a GRanges object.
        reads2 <- unlist(range(reads))
    } else {
        ## Paired-end reads. We produce a GRangesList object with 2 ranges
        ## per top-level elements.
        reads2 <- GenomicRanges:::fillGaps(reads)
    }

    sg@genes@unlistData <- .assignSubfeatureHits(reads2, ex_by_tx, ov1,
                                                 ignore.strand=TRUE,
                                                 hits.colname=hits.colname)
    sg@in_by_tx <- .assignSubfeatureHits(reads2, sg@in_by_tx, ov1,
                                         ignore.strand=TRUE,
                                         hits.colname=hits.colname)
    sg
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countReads()
###

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
    gene_id <- names(edges_by_tx)
    edges_by_tx <- sgedgesByTranscript(x, with.hits.mcols=TRUE)
    edge_data <- mcols(unlist(edge_by_tx))
    edge_data_colnames <- colnames(edge_data)
    hits_mcol_idx <- grep("\\.hits$", edge_data_colnames)
    edge_data_breakpoints <- end(PartitioningByEnd(edge_by_tx))

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

