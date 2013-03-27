\name{sgedges-methods}

\alias{sgedges-methods}

\alias{txpaths}
\alias{txpaths,SplicingGraphs-method}

\alias{UATXHcount}
\alias{UATXHcount,SplicingGraphs-method}

\alias{sgedges}
\alias{sgedges,SplicingGraphs-method}
\alias{sgedges,IntegerList-method}
\alias{sgedges,data.frame-method}

\alias{sgnodes}
\alias{sgnodes,ANY-method}
\alias{sgnodes,IntegerList-method}
\alias{sgnodes,data.frame-method}
\alias{sgnodes,DataFrame-method}

\alias{outdeg}
\alias{outdeg,ANY-method}
\alias{outdeg,DataFrame-method}

\alias{indeg}
\alias{indeg,ANY-method}
\alias{indeg,DataFrame-method}

\alias{uninformativeSSids}
\alias{uninformativeSSids,ANY-method}
\alias{uninformativeSSids,DataFrame-method}

\alias{sgedges2}


\title{
  Extract the edges (and nodes) of a splicing graph
}

\description{
  \code{sgedges} (resp. \code{sgnodes}) extracts the edges (resp. the nodes)
  of the splicing graph of a given gene from a \link{SplicingGraphs} object.
}

\usage{
sgedges(x, UATXHcount=NULL, keep.dup.edges=FALSE)
sgnodes(x)
outdeg(x)
indeg(x)

sgedges2(x)

## Related utilities:

txpaths(x, as.matrix=FALSE)
UATXHcount(x)
uninformativeSSids(x)
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object of length 1.
  }
  \item{UATXHcount}{
    TODO
  }
  \item{keep.dup.edges}{
    TODO
  }
  \item{as.matrix}{
    TODO
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

    \item \code{\link{sgedgesByTranscript}} for extracting all the splicing
          graph edges from a \link{SplicingGraphs} object.

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

sgedges(sg["geneD"])
sgnodes(sg["geneD"])
outdeg(sg["geneD"])
indeg(sg["geneD"])

txpaths(sg["geneD"])
txpaths(sg["geneD"], as.matrix=TRUE)  # splicing matrix
}
