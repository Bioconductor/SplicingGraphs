### =========================================================================
### "plotTranscripts" methods
### -------------------------------------------------------------------------

.plotTranscripts.GRangesList <- function(x, reads=NULL, from=NULL, to=NULL)
{
    ## Compute the genomic range of the transcripts.
    unlisted_x <- unlist(x, use.names=FALSE)
    strand(unlisted_x) <- "*"
    tx_range <- range(unlisted_x)
    if (length(tx_range) != 1L)
        stop("cannot plot transcripts that are on different chromosomes")

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

    if (!is.null(reads)) {
        ## Reads track.
        reads <- grglist(reads, order.as.in.query=TRUE)
        gr <- unlist(reads)
        mcols(gr)$group <- names(gr)
        reads <- relist(unname(gr), reads)
        name <- if (length(reads) == 1L) names(reads)[1L] else "reads"
        reads_track <- Gviz::AnnotationTrack(reads, name=name,
                                             fill="blue", shape="box")
        tracks <- c(tracks, list(reads_track))
    }

    if (is.null(from) || is.null(to)) {
        x_min_start <- start(tx_range)
        x_max_end <- end(tx_range)
        margin <- 0.10 * (x_max_end - x_min_start)
        if (is.null(from))
            from <- x_min_start - margin
        if (is.null(to))
            to <- x_max_end + margin
    }

    Gviz::plotTracks(tracks, from=from, to=to)
}

setGeneric("plotTranscripts", signature="x",
    function(x, reads=NULL, from=NULL, to=NULL)
        standardGeneric("plotTranscripts")
)

setMethod("plotTranscripts", "GRangesList", .plotTranscripts.GRangesList)

setMethod("plotTranscripts", "TranscriptDb",
    function(x, reads=NULL, from=NULL, to=NULL)
    {
        ex_by_tx <- exonsBy(x, by="tx", use.names=TRUE)
        plotTranscripts(ex_by_tx, reads=reads, from=from, to=to)
    }
)

setMethod("plotTranscripts", "SplicingGraphs",
    function(x)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        plotTranscripts(x[[1L]], reads=reads, from=from, to=to)
    }
)

