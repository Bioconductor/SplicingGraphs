\name{toy_data}

\alias{toy_data}
\alias{toy_genes_gff}
\alias{toy_reads_sam}
\alias{toy_reads_bam}

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
  toy_genes_txdb <- makeTranscriptDbFromGFF(toy_genes_gff())
)

## Plot the gene model:
if (interactive()) {
  library(Gviz)
  grtrack <- GeneRegionTrack(toy_genes_txdb, name="toy genes")
  gatrack <- GenomeAxisTrack()
  plotTracks(list(gatrack, grtrack))
}

## ---------------------------------------------------------------------
## Load the toy reads
## ---------------------------------------------------------------------
## The reads are single-end reads. They were aligned to the exact same
## reference genome that the above toy gene model is based on.
toy_reads_sam()
toy_reads_bam()
readGappedAlignments(toy_reads_bam(), use.names=TRUE)
}
