\name{plotTranscripts-methods}

\alias{plotTranscripts-methods}

\alias{plotTranscripts}
\alias{plotTranscripts,GRangesList-method}
\alias{plotTranscripts,TranscriptDb-method}
\alias{plotTranscripts,SplicingGraphs-method}


\title{
  Plot a set of transcripts along genomic coordinates.
}

\description{
  \code{plotTranscripts} uses the Gviz package to plot the exon structure
  of a set of transcripts along genomic coordinates.
}

\usage{
plotTranscripts(x, reads=NULL, from=NULL, to=NULL)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRangesList} object containing the genomic ranges
    of a set of exons grouped by transcript. Alternatively, \code{x} can be a
    \link[GenomicFeatures]{TranscriptDb} object, or a \link{SplicingGraphs}
    object of length 1.
  }
  \item{reads}{
    A \link[GenomicRanges]{GappedAlignments} or
    \link[GenomicRanges]{GappedAlignmentPairs} object containing single-end
    or paired-end reads.
  }
  \item{from, to}{
    TODO
  }
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \code{\link[Gviz]{plotTracks}} function in the Gviz package
          that \code{plotTranscripts} is based on.

    \item The \link[GenomicRanges]{GRangesList},
          \link[GenomicRanges]{GappedAlignments}, and
          \link[GenomicRanges]{GappedAlignmentPairs} classes
          in the GenomicRanges package.

    \item The \link[GenomicFeatures]{TranscriptDb} class in the GenomicFeatures
          package.

    \item The \link{SplicingGraphs} class.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{countReads}} for assigning reads to the edges of a
          \link{SplicingGraphs} object and summarizing them.

    \item \code{\link{toy_reads_bam}} for details about the toy data included
          in the SplicingGraphs package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. PLOT TRANSCRIPTS
## ---------------------------------------------------------------------
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

## The transcripts of a given gene can be extracted with [[. The result
## is an *unnamed* GRangesList object containing the exons grouped by
## transcript:
sg[["geneD"]]
plotTranscripts(sg[["geneD"]])  # requires the Gviz package

## The transcripts of all the genes can be extracted with unlist(). The
## result is a *named* GRangesList object containing the exons grouped
## by transcript. The names on the object are the gene ids:
ex_by_tx <- unlist(sg)
ex_by_tx
plotTranscripts(ex_by_tx)

## ---------------------------------------------------------------------
## B. PLOT TRANSCRIPTS AND READS
## ---------------------------------------------------------------------
gal <- readGappedAlignments(toy_reads_bam(), use.names=TRUE)
plotTranscripts(sg[["geneA"]], reads=gal)
plotTranscripts(ex_by_tx, reads=gal)
plotTranscripts(ex_by_tx, reads=gal, from=1, to=320)
plotTranscripts(ex_by_tx, reads=gal[21:26], from=1, to=320)
}
