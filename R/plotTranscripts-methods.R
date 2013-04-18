### =========================================================================
### "plotTranscripts" methods
### -------------------------------------------------------------------------

### Slighly faster and less memory consuming than subsetByOverlaps()
### when 'subject' contains 1 range only. Ignores the strand.
.subsetByOverlapWithSingleRange <- function(query, subject)
{
    if (!(is(query, "GAlignments") || is(query, "GAlignmentPairs")))
        stop("'query' must be a GAlignments or GAlignmentPairs object")
    if (!is(subject, "GenomicRanges"))
        stop("'subject' must be a GenomicRanges object")
    if (length(subject) != 1L)
        stop("'subject' must contain 1 range only")

    ## merge() will check that 'query' and 'subject' are based on compatible
    ## reference genomes.
    merge(seqinfo(query), seqinfo(subject))
    ## Drop all sequence levels but the one level that is used by the single
    ## range in 'subject'.
    seqlevels(subject) <- as.character(seqnames(subject))
    ## Set the same one sequence level on 'query'. Use 'force=TRUE' to remove
    ## the elements in 'query' that are not on that sequence level.
    seqlevels(query, force=TRUE) <- seqlevels(subject)

    query_ranges <- granges(query)
    strand(query_ranges) <- "*"
    cmp <- compare(query_ranges, subject)
    ## Keep reads that overlap with or are adjacent to 'subject'.
    query[-5L <= cmp & cmp <= 5L]
}

.plotTranscripts.GRangesList <- function(x, reads=NULL, from=NA, to=NA)
{
    ## Compute the genomic range of the transcripts.
    unlisted_x <- unlist(x, use.names=FALSE)
    strand(unlisted_x) <- "*"
    plot_range <- range(unlisted_x)
    if (length(plot_range) != 1L)
        stop("cannot plot transcripts that are on different chromosomes")

    ## Compute 'from' and 'to'.
    if (!(is.null(from) || isSingleNumberOrNA(from)))
        stop("'from' must be a single number, or NA, or NULL")
    if (!(is.null(to) || isSingleNumberOrNA(to)))
        stop("'to' must be a single number, or NA, or NULL")
    from_is_na <- !is.null(from) && is.na(from)
    to_is_na <- !is.null(to) && is.na(to)
    if (from_is_na || to_is_na) {
        x_min_start <- start(plot_range)
        x_max_end <- end(plot_range)
        margin <- 0.10 * (x_max_end - x_min_start)
        if (from_is_na)
            from <- x_min_start - margin
        if (to_is_na)
            to <- x_max_end + margin
        from_str <- ifelse(is.null(from), "NULL", from)
        to_str <- ifelse(is.null(to), "NULL", to)
        message("using: from=", from_str, ", to=", to_str)
    }

    ## Genome axis.
    tracks <- list(Gviz::GenomeAxisTrack())

    ## Transcript tracks (we create 1 track per transcript).
    track_names <- mcols(x)$tx_id
    if (is.null(track_names))
        track_names <- names(x)
    tx_tracks <- lapply(seq_along(x),
                        function(i) {
                          tx <- x[[i]]
                          Gviz::AnnotationTrack(tx, name=track_names[i],
                                                fill="orange", shape="box")
                        })
    tracks <- c(tracks, tx_tracks)

    ## Reads track.
    if (!is.null(reads)) {
        if (length(reads) != 0L && !(is.null(from) && is.null(to))) {
            if (is.null(from)) {
                start(plot_range) <- min(start(reads))
            } else {
                start(plot_range) <- from
            }
            if (is.null(to)) {
                end(plot_range) <- max(end(reads))
            } else {
                end(plot_range) <- to
            }
            reads <- .subsetByOverlapWithSingleRange(reads, plot_range)
        }
        reads <- grglist(reads, order.as.in.query=TRUE)

        ## Set strand to *.
        strand(reads@unlistData) <- "*"

        ## Set group id (required by the Gviz package, why?)
        mcols(reads@unlistData)$group <- rep.int(seq_along(reads), 
                                                 elementLengths(reads))
        name <- if (length(reads) == 1L) names(reads)[1L] else "reads"
        reads_track <- Gviz::AnnotationTrack(reads, name=name,
                                             fill="blue", shape="box")
        tracks <- c(tracks, list(reads_track))
    }

    Gviz::plotTracks(tracks, from=from, to=to)
}

setGeneric("plotTranscripts", signature="x",
    function(x, reads=NULL, from=NA, to=NA)
        standardGeneric("plotTranscripts")
)

setMethod("plotTranscripts", "GRangesList", .plotTranscripts.GRangesList)

setMethod("plotTranscripts", "TranscriptDb",
    function(x, reads=NULL, from=NA, to=NA)
    {
        ex_by_tx <- exonsBy(x, by="tx", use.names=TRUE)
        plotTranscripts(ex_by_tx, reads=reads, from=from, to=to)
    }
)

setMethod("plotTranscripts", "SplicingGraphs",
    function(x, reads=NULL, from=NA, to=NA)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        plotTranscripts(x[[1L]], reads=reads, from=from, to=to)
    }
)

