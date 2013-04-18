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
plotTranscripts(x, reads=NULL, from=NA, to=NA)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRangesList} object containing the genomic ranges
    of a set of exons grouped by transcript. Alternatively, \code{x} can be a
    \link[GenomicFeatures]{TranscriptDb} object, or a \link{SplicingGraphs}
    object of length 1.
  }
  \item{reads}{
    A \link[GenomicRanges]{GAlignments} or
    \link[GenomicRanges]{GAlignmentPairs} object containing single-end
    or paired-end reads.
  }
  \item{from, to}{
    Single numeric values, giving the range of genomic coordinates to plot
    the tracks in. By default (i.e. \code{from=NA} and \code{to=NA}), the
    plot covers the range spanned by the transcripts.
    If \code{from=NULL} and \code{to=NULL}, then the plot covers the
    range spanned by the transcripts and the reads.
  }
}

\author{
  H. Pages
}

\seealso{
  This man page is part of the SplicingGraphs package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.

  Other topics related to this man page and documented in other packages:
  \itemize{
    \item The \code{\link[Gviz]{plotTracks}} function in the Gviz package
          that \code{plotTranscripts} is based on.

    \item The \link[GenomicRanges]{GRangesList},
          \link[GenomicRanges]{GAlignments}, and
          \link[GenomicRanges]{GAlignmentPairs} classes
          in the GenomicRanges package.

    \item The \link[GenomicFeatures]{TranscriptDb} class in the GenomicFeatures
          package.
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
gal <- readGAlignments(toy_reads_bam(), use.names=TRUE)
plotTranscripts(sg[["geneA"]], reads=gal)
plotTranscripts(ex_by_tx, reads=gal)
plotTranscripts(ex_by_tx, reads=gal, from=1, to=320)
plotTranscripts(ex_by_tx, reads=gal[21:26], from=1, to=320)
}
