\name{sgedges-methods}

\alias{sgedges-methods}

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

## Related utility:
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

    \item \code{\link{sgedgesByTranscript}} for extracting all the introns
          or splicing graph edges from a \link{SplicingGraphs} object.

    \item \code{\link{txpath}} for extracting the transcript paths of a
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

sgedges(sg["geneD"])
sgnodes(sg["geneD"])
outdeg(sg["geneD"])
indeg(sg["geneD"])

## Sanity checks:
check_way1_vs_way2 <- function(res1, res2)
{
    edges1 <- res1[res1$ex_or_in != "", ]  # remove artificial edges
    edges2 <- mcols(unlist(res2, use.names=FALSE))
    edges2 <- edges2[ , -1]  # remove "sgedge_id" col (global edge id)
    stopifnot(identical(edges1, edges2))
}

for (i in seq_along(sg)) {
    sgi <- sg[i]
    check_way1_vs_way2(
        ## After removal of the artificial edges, the edges returned
        ## by 'sgedges( , keep.dup.edges=TRUE)' should be the same as
        ## those returned by 'sgedgesByTranscript()' on a SplicingGraphs
        ## object of length 1.
        DataFrame(sgedges(sgi, keep.dup.edges=TRUE)),
        sgedgesByTranscript(sgi))
    check_way1_vs_way2(
        ## After removal of the artificial edges, the edges returned
        ## by 'sgedges()' should be the same as those returned
        ## by 'sgedgesByGene()' on a SplicingGraphs object of length 1.
        sgedges(sgi),
        sgedgesByGene(sgi))
}
}
