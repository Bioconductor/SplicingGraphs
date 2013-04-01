### =========================================================================
### Functions for assigning reads to the edges of a SplicingGraphs object and
### for summarizing them
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .assignSubfeatureHits()
###
### This is the workhorse behind assignReads().
###

### 'query': a named GRangesList object containing gapped reads.
###     Alternatively it can be the GRanges object obtained by extracting
###     a single range per read e.g. the spanning ranges obtained with
###     'unlist(range(query))'. In any case, it must have the length and
###     names of the original GRangesList object.
### 'subject': a GRangesList object containing the ranges of some subfeature
###     (e.g. exonic or intronic ranges) grouped by their parent feature
###     (e.g. transcripts).
### 'hits': a Hits object compatible with 'query' and 'subject'.
.check_assignSubfeatureHits_args <- function(query, subject, hits,
                                             ignore.strand, hits.colname)
{
    if (!(is(query, "GRangesList") || is(query, "GRanges")))
        stop("'query' must be a GRangesList or GRanges object")
    query_names <- names(query)
    if (is.null(query_names))
        stop("'query' must have names")
    if (anyDuplicated(query_names))
        stop("'query' has duplicated names")
    if (!is(subject, "GRangesList"))
        stop("'subject' must be a GRangesList object")
    if (!is(hits, "Hits"))
        stop("'hits' must be a Hits object")
    if (queryLength(hits) != length(query)
     || subjectLength(hits) !=  length(subject))
        stop("'hits' is not compatible with 'query' and 'subject'")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    if (!isSingleString(hits.colname))
        stop("'hits.colname' must be a single string")
}

### Returns 'subject' with 1 additional inner metadata col "hits" containing
### the hits assigned to each subrange.
.assignSubfeatureHits <- function(query, subject, hits, ignore.strand=FALSE,
                                  hits.colname="hits")
{
    .check_assignSubfeatureHits_args(query, subject, hits,
                                     ignore.strand, hits.colname)
    query_names <- names(query)
    subject_unlistData <- subject@unlistData
    subhits <- findOverlaps(query, subject_unlistData,
                            ignore.strand=ignore.strand)
    subhits_q <- queryHits(subhits)
    subhits_s <- togroup(subject@partitioning, subjectHits(subhits))
    m <- IRanges:::matchIntegerPairs(subhits_q, subhits_s,
                                     queryHits(hits), subjectHits(hits))
    subhits <- subhits[!is.na(m)]
    hit_per_subfeature <- splitAsList(query_names[queryHits(subhits)],
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

### FIXME: It's questionable whether this does the right thing on paired-end
### reads. I guess not...
assignReads <- function(sg, reads, sample.name=NA)
{
    if (!is(sg, "SplicingGraphs"))
        stop("'sg' must be a SplicingGraphs object")
    if (is(reads, "GappedAlignments")) {
        reads <- grglist(reads, order.as.in.query=TRUE)
    } else if (!is(reads, "GRangesList")) {
        stop("'reads' must be a GRangesList object")
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
    ov0_is_comp <- isCompatibleWithSplicing(ovenc0)
    ov1 <- ov0[ov0_is_comp]
    reads2 <- unlist(range(reads))
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

### Return a DataFrame with 1 row per unique splicing graph edge and 1 column
### per sample.
countReads <- function(sg)
{
    if (!is(sg, "SplicingGraphs"))
        stop("'sg' must be a SplicingGraphs object")
    edges_by_gene <- sgedgesByGene(sg, with.hits.mcols=TRUE)
    edges0 <- unlist(edges_by_gene, use.names=FALSE)
    edges0_mcols <- mcols(edges0)
    edges0_mcolnames <- colnames(edges0_mcols)
    hits_idx <- grep("\\.hits$", edges0_mcolnames)
    ans <- endoapply(edges0_mcols[hits_idx], elementLengths)
    colnames(ans) <- sub("\\.hits$", "", colnames(ans))
    left_cols <- edges0_mcols[ , c("sgedge_id", "ex_or_in")]
    cbind(left_cols, ans)
}

