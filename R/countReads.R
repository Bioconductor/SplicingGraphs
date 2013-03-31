### =========================================================================
### Functions for counting compatible hits per transcript, and for assigning
### compatible hits per exon or per intron
### -------------------------------------------------------------------------


### 'query': a named GRangesList object containing gapped reads.
### 'subject': a GRangesList object containing some subfeature (e.g. exons
###     or introns) grouped by their parent feature (e.g. transcripts).
### 'hits': a Hits object compatible with 'query' and 'subject'.
### Returns 'subject' with additional inner metadata col "hits"
### (CharacterList) reporting the hits for each subfeature.
### TODO: Current implementation is messy and inefficient. There must be
### a better way...
.assignSubfeatureHits <- function(query, subject, hits, ignore.strand=FALSE,
                                  hits.colname="hits")
{
    if (!is(query, "GRangesList"))
        stop("'query' must be a GRangesList object")
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
    unlisted_subject <- subject@unlistData
    #if (hits.colname %in% colnames(unlisted_subject))
    #    stop("'unlisted(subject)' already has metadata column ", hits.colname)

    subject_eltlens <- elementLengths(subject)  # nb of subfeatures per subject
    s_hits <- subjectHits(hits)
    tx1 <- subject[s_hits]
    tx1_eltlens <- subject_eltlens[s_hits]
    ex11 <- unlist(tx1, use.names=FALSE)

    q_hits <- queryHits(hits)
    q_hits11 <- rep.int(q_hits, tx1_eltlens)
    gr11 <- unlist(range(query), use.names=FALSE)[q_hits11]

    if (ignore.strand)
        strand(gr11) <- strand(ex11) <- "*"

    cmp11 <- compare(gr11, ex11)
    is_hit11 <- -4L <= cmp11 & cmp11 <= 4L
    ex11_hit <- ifelse(is_hit11, query_names[q_hits11], NA_character_)
    mcols(tx1@unlistData)$hit <- ex11_hit

    unq_s_hits <- unique(s_hits)
    subfeature_hits <- lapply(unq_s_hits,
                         function(i) {
                           nsubfeatures <- subject_eltlens[i]
                           hits <- splitAsList(
                                     mcols(tx1[s_hits == i]@unlistData)$hit,
                                     seq_len(nsubfeatures))
                           if (length(hits) != 0L)
                               hits <- hits[!is.na(hits)]
                           hits
                         })

    mcols(subject@unlistData)[[hits.colname]] <- CharacterList(character(0))
    mcols(subject[unq_s_hits]@unlistData)[[hits.colname]] <-
                                                 do.call(c, subfeature_hits)
    subject
}

### FIXME: It's questionable whether this does the right thing on paired-end
### reads. I guess not...
assignReads <- function(sg, reads, sample.name=NA)
{
    if (!is(sg, "SplicingGraphs"))
        stop("'sg' must be a SplicingGraphs object")
    if (!is(reads, "GRangesList"))
        stop("'reads' must be a GRangesList object")
    if (!isSingleStringOrNA(sample.name))
        stop("'sample.name' must be a single string or NA")
    if (is.na(sample.name)) {
        hits.colname <- "hits"
    } else {
        hits.colname <- paste0(sample.name, ".hits")
    }

    unlisted_sg <- unlist(sg)
    ov0 <- findOverlaps(reads, unlisted_sg, ignore.strand=TRUE)
    ovenc0 <- encodeOverlaps(reads, unlisted_sg, hits=ov0,
                             flip.query.if.wrong.strand=TRUE)
    ov0_is_comp <- isCompatibleWithSplicing(ovenc0)
    ov1 <- ov0[ov0_is_comp]
    sg@genes@unlistData <- unname(.assignSubfeatureHits(reads, unlisted_sg,
                                                    ov1, ignore.strand=TRUE,
                                                    hits.colname=hits.colname))
    sg@in_by_tx <- .assignSubfeatureHits(reads, sg@in_by_tx, ov1,
                                         ignore.strand=TRUE,
                                         hits.colname=hits.colname)
    sg
}

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

