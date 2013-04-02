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
plotTranscripts(x)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRangesList} object containing the genomic ranges
    of a set of exons grouped by transcript. Alternatively, \code{x} can be a
    \link[GenomicFeatures]{TranscriptDb} object, or a \link{SplicingGraphs}
    object of length 1.
  }
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \code{\link[Gviz]{plotTracks}} function in the Gviz package
          that \code{plotTranscripts} is based on.

    \item The \link[GenomicRanges]{GRangesList} class in the GenomicRanges
          package.

    \item The \link[GenomicFeatures]{TranscriptDb} class in the GenomicFeatures
          package.

    \item The \link{SplicingGraphs} class.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.
  }
}

\examples{
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
}
