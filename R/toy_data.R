### =========================================================================
### Little helpers for quick access to the toy data
### -------------------------------------------------------------------------

toy_genes_gff <- function()
{
    system.file("extdata", "toy_genes.gff3",
                package="SplicingGraphs", mustWork=TRUE)
}

toy_reads_sam <- function()
{
    system.file("extdata", "toy_reads.sam",
                package="SplicingGraphs", mustWork=TRUE)
}

.toy_reads_bam_cache <- new.env(parent=emptyenv())

toy_reads_bam <- function()
{
    toy_reads_bam <- try(get("toy_reads_bam", envir=.toy_reads_bam_cache,
                             inherits=FALSE), silent=TRUE)
    if (!is(toy_reads_bam, "try-error"))
        return(toy_reads_bam)
    toy_reads_sam <- toy_reads_sam()
    destination <- tempfile()
    toy_reads_bam <- asBam(toy_reads_sam, destination)
    ## Should never happen.
    if (toy_reads_bam != paste0(destination, ".bam"))
        stop("asBam() returned an unexpected path")
    assign("toy_reads_bam", toy_reads_bam, envir=.toy_reads_bam_cache)
    toy_reads_bam
}

### Displaying bug in Gviz (reported to Florian on Feb. 15, 2013). The first
### range has width 5 and therefore should cover 5 letters (including the 1st
### T). The last range has width 1 (not 0) and should be visible.
if (FALSE) {
  library(Biostrings)
  seq <- DNAStringSet(c(chrX="ACCGACTTCA"))
  library(GenomicRanges)
  gr <- GRanges("chrX", IRanges(3, width=5:1))
  library(Gviz) 
  plotTracks(list(GenomeAxisTrack(), SequenceTrack(seq), AnnotationTrack(gr)),
             from=1, to=9)
}

plotToyReads <- function(gal, txdb, from=NULL, to=NULL)
{
    ax_track <- GenomeAxisTrack()
    txdb_track <- GeneRegionTrack(txdb, name="genes")
    grl <- grglist(gal)
    gr <- unlist(grl)
    mcols(gr)$group <- names(gr)
    grl <- relist(unname(gr), grl)
    name <- if (length(grl) == 1L) names(grl)[1L] else "reads"
    gal_track <- AnnotationTrack(grl, name=name, fill="blue", shape="box")
    if (is.null(from) || is.null(to)) {
        gal_min_start <- min(start(gal))
        gal_max_end <- max(end(gal))
        margin <- 0.10 * (gal_max_end - gal_min_start)
        if (is.null(from))
            from <- gal_min_start - margin
        if (is.null(to))
            to <- gal_max_end + margin
    }
    plotTracks(list(ax_track, gal_track, txdb_track), from=from, to=to)
}

