\name{toy_data}

\alias{toy_data}
\alias{toy_genes_gff}
\alias{toy_reads_sam}
\alias{toy_reads_bam}
\alias{toy_overlaps}

\title{
  Little helpers for quick access to the toy data included in the
  \pkg{SplicingGraphs} package
}

\description{
  TODO
}

\usage{
toy_genes_gff()
toy_reads_sam()
toy_reads_bam()
toy_overlaps()
}

\author{
  H. Pages
}

\seealso{
  This man page is part of the \pkg{SplicingGraphs} package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.

  Other topics related to this man page and documented in other packages:
  \itemize{
    \item The \link[GenomicRanges]{GRangesList} class defined in the
          \pkg{GenomicRanges} package.

    \item The \link[GenomicAlignments]{GAlignments} and
          \link[GenomicAlignments]{GAlignmentPairs} classes defined
          in the \pkg{GenomicAlignments} package.

    \item The \code{\link[GenomicFeatures]{makeTxDbFromGFF}}
          function and the \link[GenomicFeatures]{TxDb} class
          defined in the \pkg{GenomicFeatures} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. LOAD THE TOY GENE MODEL AS A TxDb OBJECT AND PLOT IT
## ---------------------------------------------------------------------
toy_genes_gff()

## Note that you can display the content of the file with:
cat(readLines(toy_genes_gff()), sep="\n")

library(GenomicFeatures)
suppressWarnings(
  txdb <- makeTxDbFromGFF(toy_genes_gff())
)

## Plot all the transcripts in the gene model:
plotTranscripts(txdb)

## ---------------------------------------------------------------------
## B. LOAD THE TOY READS AS A GAlignments OBJECT AND PLOT THEM
## ---------------------------------------------------------------------
## The reads are single-end reads. They are assumed to come from an
## RNA-seq experiment and to have been aligned to the exact same
## reference genome that the above toy gene model is based on.
toy_reads_sam()
toy_reads_bam()
gal <- readGAlignments(toy_reads_bam(), use.names=TRUE)

plotTranscripts(txdb, reads=gal)
plotTranscripts(txdb, reads=gal, from=1, to=320)

## ---------------------------------------------------------------------
## C. FIND THE OVERLAPS BETWEEN THE TOY READS AND THE TOY GENE MODEL
## ---------------------------------------------------------------------
grl <- grglist(gal, order.as.in.query=TRUE)
ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)

## Most of the times the RNA-seq protocol is unstranded so the strand
## reported in the BAM file for each alignment is meaningless. Thus we
## should call findOverlaps() with 'ignore.strand=TRUE':
ov0 <- findOverlaps(grl, ex_by_tx, ignore.strand=TRUE)

## Put the overlaps in a data.frame to make it easier to read:
df0 <- data.frame(QNAME=names(grl)[queryHits(ov0)],
                  tx_id=names(ex_by_tx)[subjectHits(ov0)],
                  stringsAsFactors=FALSE)
head(df0)

## These overlaps have been manually checked and included in the
## SplicingGraphs package. They can be loaded with the toy_overlaps()
## helper:
toy_ov <- toy_overlaps()
head(toy_ov)
stopifnot(identical(df0, toy_ov[ , 1:2]))

## ---------------------------------------------------------------------
## D. DETECT THE OVERLAPS THAT ARE COMPATIBLE WITH THE GENE MODEL
## ---------------------------------------------------------------------
## First we encode the overlaps:
ovenc0 <- encodeOverlaps(grl, ex_by_tx, hits=ov0,
                         flip.query.if.wrong.strand=TRUE)
ovenc0

## Each encoding tells us whether the corresponding overlap is
## compatible or not with the gene model:
ov0_is_comp <- isCompatibleWithSplicing(ovenc0)
head(ov0_is_comp)

## Overlap compatibility has also been manually checked and included in
## the table returned by toy_overlaps():
stopifnot(identical(ov0_is_comp, toy_ov[ , 3]))
}
