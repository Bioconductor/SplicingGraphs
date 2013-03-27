\name{sgedgesByTranscript-methods}

\alias{sgedgesByTranscript-methods}

\alias{sgedgesByTranscript}
\alias{sgedgesByTranscript,SplicingGraphs-method}


\title{
  Extract all the splicing graph edges from a SplicingGraphs object
}

\description{
  \code{sgedgesByTranscript} extracts the splicing graph edges of all
  the genes from a \link{SplicingGraphs} object. Return them in a
  \link[GenomicRanges]{GRangesList} object where they are grouped by
  transcript.
}

\usage{
sgedgesByTranscript(x)
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object.
  }
}

\details{
  TODO
}

\value{
  TODO
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \link{SplicingGraphs} class.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{bubbles}} for extracting the bubbles of a splicing graph.
  }
}

\examples{
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

ex_by_tx <- unlist(sg)
in_by_tx <- intronsByTranscript(sg)

## 'sgedgesByTranscript(sg)' combines 'ex_by_tx' and 'in_by_tx' in
## a single GRangesList object.
sgedges_by_tx <- sgedgesByTranscript(sg)

## A sanity check:
stopifnot(identical(elementLengths(sgedges_by_tx),
                    elementLengths(ex_by_tx) + elementLengths(in_by_tx)))
}
