### =========================================================================
### Functions for counting compatible hits per transcript, and for assigning
### compatible hits per exon or per intron
### -------------------------------------------------------------------------

setMethod("findOverlaps", c("GRangesList", "SplicingGraphs"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=ignore.strand)
    {
        findOverlaps(query, subject@tx,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

setMethod("encodeOverlaps", c("GRangesList", "SplicingGraphs"),
    function(query, subject, hits=NULL, flip.query.if.wrong.strand=FALSE)
    {
        encodeOverlaps(query, subject@tx,
                       hits=hits,
                       flip.query.if.wrong.strand=flip.query.if.wrong.strand)
    }
)

### 'query': a named GRangesList object containing gapped reads.
### 'subject': a GRangesList object containing some subfeature (e.g. exons
###     or introns) grouped by their parent feature (e.g. transcripts).
### 'hits': a Hits object compatible with 'query' and 'subject'.
### Returns 'subject' with additional inner metadata col "hits"
### (CharacterList) reporting the hits for each subfeature.
### TODO: Current implementation is messy and inefficient. There must be
### a better way...
assignSubfeatureHits <- function(query, subject, hits, ignore.strand=FALSE)
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

    mcols(subject@unlistData)$hits <- CharacterList(character(0))
    mcols(subject[unq_s_hits]@unlistData)$hits <- do.call(c, subfeature_hits)
    subject
}

