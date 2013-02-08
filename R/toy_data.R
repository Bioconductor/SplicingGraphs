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

