\name{sgedgesByTranscript-methods}

\alias{sgedgesByTranscript-methods}

\alias{intronsByTranscript,SplicingGraphs-method}

\alias{sgedgesByTranscript}
\alias{sgedgesByTranscript,SplicingGraphs-method}


\title{
  Extract all the introns or splicing graph edges from a SplicingGraphs object
}

\description{
  \code{intronsByTranscript} extracts all the introns of all
  the genes from a \link{SplicingGraphs} object.

  \code{sgedgesByTranscript} extracts all the splicing graph edges of all
  the genes from a \link{SplicingGraphs} object.

  Both functions return a named \link[GenomicRanges]{GRangesList} object
  where the items (i.e. introns for the former, splicing graph edges for
  the latter) are grouped by transcript. The names of this object are the
  gene ids.
}

\usage{
\S4method{intronsByTranscript}{SplicingGraphs}(x)

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

    \item \code{\link{txpath}} for extracting the transcript paths of a
          splicing graph.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{bubbles}} for computing the bubbles of a splicing graph.

    \item \code{\link{countReads}} for assigning reads to a
          \link{SplicingGraphs} object and counting them.
  }
}

\examples{
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

in_by_tx <- intronsByTranscript(sg)
in_by_tx

## A sanity check:
ex_by_tx <- unlist(sg)
stopifnot(identical(elementLengths(in_by_tx) + 1L,
                    elementLengths(ex_by_tx)))

## 'sgedgesByTranscript(sg)' combines 'ex_by_tx' and 'in_by_tx' in
## a single GRangesList object.
sgedges_by_tx <- sgedgesByTranscript(sg)

## A sanity check:
stopifnot(identical(elementLengths(sgedges_by_tx),
                    elementLengths(ex_by_tx) + elementLengths(in_by_tx)))
}
