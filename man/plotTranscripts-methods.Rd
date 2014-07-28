\name{plotTranscripts-methods}

\alias{plotTranscripts-methods}

\alias{plotTranscripts}
\alias{plotTranscripts,GRangesList-method}
\alias{plotTranscripts,TxDb-method}
\alias{plotTranscripts,SplicingGraphs-method}


\title{
  Plot a set of transcripts along genomic coordinates.
}

\description{
  \code{plotTranscripts} uses the \pkg{Gviz} package to plot the exon
  structure of a set of transcripts along genomic coordinates.
}

\usage{
plotTranscripts(x, reads=NULL, from=NA, to=NA, max.plot.reads=200)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRangesList} object containing the genomic ranges
    of a set of exons grouped by transcript. Alternatively, \code{x} can be a
    \link[GenomicFeatures]{TxDb} object, or a \link{SplicingGraphs}
    object of length 1.
  }
  \item{reads}{
    A \link[GenomicAlignments]{GAlignments} or
    \link[GenomicAlignments]{GAlignmentPairs} object containing single-end
    or paired-end reads.
  }
  \item{from, to}{
    Single numeric values, giving the range of genomic coordinates to plot
    the tracks in. By default (i.e. \code{from=NA} and \code{to=NA}), the
    plot covers the range spanned by the transcripts.
    If \code{from=NULL} and \code{to=NULL}, then the plot covers the
    range spanned by the transcripts and the reads.
  }
  \item{max.plot.reads}{
    The maximum number of reads that will be plotted. When the number of
    reads that fall in the region being plotted is very large, plotting them
    all would take a long time and result in a plot that is not very useful.
    If that number is greater than \code{max.plot.reads}, then only
    \code{max.plot.reads} randomly chosen reads are plotted.
  }
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
    \item \code{plotTranscripts} is based on the \code{\link[Gviz]{plotTracks}}
          function defined in the \pkg{Gviz} package.

    \item The \link[GenomicRanges]{GRangesList} class defined in the
          \pkg{GenomicRanges} package.

    \item The \link[GenomicAlignments]{GAlignments} and
          \link[GenomicAlignments]{GAlignmentPairs} classes
          defined in the \pkg{GenomicAlignments} package.

    \item The \link[GenomicFeatures]{TxDb} class defined in the
          \pkg{GenomicFeatures} package.
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
