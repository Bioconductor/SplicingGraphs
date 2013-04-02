### =========================================================================
### "plotTranscripts" methods
### -------------------------------------------------------------------------


setGeneric("plotTranscripts", function(x) standardGeneric("plotTranscripts"))

setMethod("plotTranscripts", "GRangesList",
    function(x)
    {
        track_names <- mcols(x)$tx_id
        if (is.null(track_names))
            track_names <- names(x)
        ### We create 1 track per transcript.
        tx_tracks <- lapply(seq_along(x),
                            function(i) {
                              tx <- x[[i]]
                              Gviz::AnnotationTrack(tx, name=track_names[i],
                                                    fill="orange", shape="box")
                            })
        ax_track <- Gviz::GenomeAxisTrack()
        Gviz::plotTracks(c(list(ax_track), tx_tracks))
    }
)

setMethod("plotTranscripts", "TranscriptDb",
    function(x)
    {
        ex_by_tx <- exonsBy(x, by="tx", use.names=TRUE)
        plotTranscripts(ex_by_tx)
    }
)

setMethod("plotTranscripts", "SplicingGraphs",
    function(x)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        plotTranscripts(x[[1L]])
    }
)

