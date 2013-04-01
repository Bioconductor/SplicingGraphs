\name{txpath-methods}

\alias{txpath-methods}

\alias{txpath}
\alias{txpath,GRangesList-method}
\alias{txpath,SplicingGraphs-method}

\alias{UATXHcount}
\alias{UATXHcount,SplicingGraphs-method}


\title{
  Extract the transcript paths of a splicing graph
}

\description{
  \code{txpath} extracts the transcript paths of the splicing graph
  of a given gene from a \link{SplicingGraphs} object.
}

\usage{
txpath(x, as.matrix=FALSE)

## Related utility:
UATXHcount(x)
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object of length 1
    or a \link[GenomicRanges]{GRangesList} object.
  }
  \item{as.matrix}{
    TODO
  }
}

\details{
  TODO
}

\value{
  A named list-like object with one list element per transcript in the gene.
  Each list element is an integer vector that describes the \emph{path}
  of the transcript i.e. the \emph{Splicing Site ids} that it goes thru.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \link{SplicingGraphs} class.

    \item \code{\link{sgedgesByGene}} for extracting all the edges and
          ranges from a \link{SplicingGraphs} object.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{bubbles}} for computing the bubbles of a splicing graph.

    \item \code{\link{countReads}} for assigning reads to the edges of a
          \link{SplicingGraphs} object and summarizing them.
  }
}

\examples{
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

## Note that the list elements in the returned IntegerList object
## always consist of an even number of Splicing Site ids in ascending
## order.
txpath(sg["geneB"])
txpath(sg["geneD"])
strand(sg)

txpath(sg["geneD"], as.matrix=TRUE)  # splicing matrix
}
