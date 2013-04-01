\name{toy_data}

\alias{toy_data}
\alias{toy_genes_gff}
\alias{toy_reads_sam}
\alias{toy_reads_bam}
\alias{toy_overlaps}
\alias{plotToyReads}

\title{
  Little helpers for quick access to (and plotting of) the toy data
  included in the SplicingGraphs package
}

\description{
  TODO
}

\usage{
toy_genes_gff()
toy_reads_sam()
toy_reads_bam()
toy_overlaps()
plotToyReads(gal, txdb, from=NULL, to=NULL)
}

\arguments{
  \item{gal}{
    TODO
  }
  \item{txdb}{
    TODO
  }
  \item{from}{
    TODO
  }
  \item{to}{
    TODO
  }
}

\author{
  H. Pages
}

\seealso{
  The toy data described here is used in the examples of the following
  man pages:
  \itemize{
    \item \code{\link{SplicingGraphs}}: Tools for creating and
          plotting splicing graphs from transcript annotations.

    \item \code{\link{countReads}} for assigning reads to the edges of a
          \link{SplicingGraphs} object and summarizing them.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. LOAD THE TOY GENE MODEL AS A TranscriptDb OBJECT AND PLOT IT
## ---------------------------------------------------------------------
toy_genes_gff()

## Note that you can display the content of the file with:
cat(readLines(toy_genes_gff()), sep="\n")

library(GenomicFeatures)
suppressWarnings(
  txdb <- makeTranscriptDbFromGFF(toy_genes_gff())
)

## Plot all the transcripts in the gene model:
if (interactive()) {
  plotTranscripts(txdb)
}

## ---------------------------------------------------------------------
## B. LOAD THE TOY READS AS A GappedAlignments OBJECT AND PLOT THEM
## ---------------------------------------------------------------------
## The reads are single-end reads. They are assumed to come from an
## RNA-seq experiment and to have been aligned to the exact same
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

## ---------------------------------------------------------------------
## C. FIND THE OVERLAPS BETWEEN THE TOY READS AND THE TOY GENE MODEL
## ---------------------------------------------------------------------
grl <- grglist(gal, order.as.in.query=TRUE)
ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)

## Most of the times the RNA-seq protocol was unstranded so the strand
## reported in the BAM file for each alignment is meaningless. Thus we
## call findOverlaps() with 'ignore.strand=TRUE':
ov <- findOverlaps(grl, ex_by_tx, ignore.strand=TRUE)

## Put the overlaps in a data.frame to make it easier to read:
ov2 <- data.frame(QNAME=names(grl)[queryHits(ov)],
                  tx_id=names(ex_by_tx)[subjectHits(ov)],
                  stringsAsFactors=FALSE)
head(ov2)

## These overlaps have been manually checked and included in the
## SplicingGraphs package. They can be loaded with the toy_overlaps()
## helper:
toy_ov <- toy_overlaps()
head(toy_ov)
stopifnot(identical(ov2, toy_ov[ , 1:2]))

## ---------------------------------------------------------------------
## D. DETECT THE OVERLAPS THAT ARE COMPATIBLE WITH THE GENE MODEL
## ---------------------------------------------------------------------
## First we encode the overlaps:
ovenc <- encodeOverlaps(grl, ex_by_tx, hits=ov,
                        flip.query.if.wrong.strand=TRUE)
ovenc

## Each encoding tells us whether the corresponding overlap is
## compatible or not with the gene model:
is_comp <- isCompatibleWithSplicing(ovenc)
head(is_comp)

## Overlap compatibility has also been manually checked and included in
## the table returned by toy_overlaps():
stopifnot(identical(is_comp, toy_ov[ , 3]))
}
