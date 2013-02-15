\name{toy_data}

\alias{toy_data}
\alias{toy_genes_gff}
\alias{toy_reads_sam}
\alias{toy_reads_bam}
\alias{plotToyReads}

\title{
  Little helpers for quick access to the toy data included in the
  SplicingGraphs package
}

\description{
  TODO
}

\usage{
toy_genes_gff()
toy_reads_sam()
toy_reads_bam()
}

\author{
  H. Pages
}

\seealso{
  The toy data described here is used in the examples of the following
  man pages:
  \itemize{
    \item \code{\link{splicingGraphs}}: Tools for creating and
          plotting splicing graphs from transcript annotations.

    \item \link{countReads}: Tools for assigning read counts to splicing
          graphs.
  }
}

\examples{
## ---------------------------------------------------------------------
## Load the toy gene model as a TranscriptDb object and plot it
## ---------------------------------------------------------------------
toy_genes_gff()

## Note that you can display the content of the file with:
cat(readLines(toy_genes_gff()), sep="\n")

library(GenomicFeatures)
suppressWarnings(
  txdb <- makeTranscriptDbFromGFF(toy_genes_gff())
)

## Plot the gene model:
if (interactive()) {
  library(Gviz)
  ax_track <- GenomeAxisTrack()
  tg_track <- GeneRegionTrack(txdb, name="toy genes")
  plotTracks(list(ax_track, tg_track))
}

## ---------------------------------------------------------------------
## Load the toy reads
## ---------------------------------------------------------------------
## The reads are single-end reads. They were aligned to the exact same
## reference genome that the above toy gene model is based on.
toy_reads_sam()
toy_reads_bam()
gal <- readGappedAlignments(toy_reads_bam(), use.names=TRUE)

## Plot the reads:
if (interactive()) {
  ## All of them:
  plotToyReads(gal, txdb)
  ## The first 11 only:
  plotToyReads(head(gal, n=11), txdb)
  ## The rest of them:
  plotToyReads(tail(gal, n=-11), txdb)
  ## Reads 21 to 26:
  plotToyReads(gal[21:26], txdb)
  ## Read 26 only:
  plotToyReads(gal[26], txdb)
}
