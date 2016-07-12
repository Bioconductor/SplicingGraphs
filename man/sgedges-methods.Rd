\name{sgedges-methods}

\alias{sgedges-methods}

\alias{sgedges}
\alias{sgedges,SplicingGraphs-method}

\alias{sgnodes}
\alias{sgnodes,SplicingGraphs-method}
\alias{sgnodes,IntegerList-method}
\alias{sgnodes,data.frame-method}
\alias{sgnodes,DataFrame-method}

\alias{outdeg}
\alias{outdeg,ANY-method}
\alias{outdeg,DataFrame-method}

\alias{indeg}
\alias{indeg,ANY-method}
\alias{indeg,DataFrame-method}


\title{
  Extract the edges (and nodes) of a splicing graph
}

\description{
  \code{sgedges} (resp. \code{sgnodes}) extracts the edges (resp. the nodes)
  of the splicing graph of a given gene from a \link{SplicingGraphs} object.
}

\usage{
sgedges(x, txweight=NULL, keep.dup.edges=FALSE)
sgnodes(x)
outdeg(x)
indeg(x)
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object of length 1.
  }
  \item{txweight}{
    TODO
  }
  \item{keep.dup.edges}{
    If \code{FALSE} (the default), then edges with the same \emph{global
    edge id} are merged into a single row. Use \code{keep.dup.edges=TRUE}
    if this merging should not be performed.
  }
}

\details{
  TODO
}

\value{
  TODO
}

\author{
  H. Pagès
}

\seealso{
  This man page is part of the \pkg{SplicingGraphs} package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.
}

\examples{
## ---------------------------------------------------------------------
## 1. Make SplicingGraphs object 'sg' from toy gene model (see
##    '?SplicingGraphs')
## ---------------------------------------------------------------------
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

## ---------------------------------------------------------------------
## 2. Basic usage
## ---------------------------------------------------------------------
sgedges(sg["geneD"])
sgnodes(sg["geneD"])
outdeg(sg["geneD"])
indeg(sg["geneD"])

## ---------------------------------------------------------------------
## 3. Sanity checks
## ---------------------------------------------------------------------
check_way1_vs_way2 <- function(res1, res2)
{
    edges1 <- res1[res1$ex_or_in != "", ]  # remove artificial edges
    edges2 <- mcols(unlist(res2, use.names=FALSE))
    stopifnot(identical(edges1, edges2))
}

for (i in seq_along(sg)) {
    sgi <- sg[i]
    ## After removal of the artificial edges, the edges returned
    ## by 'sgedges()' should be the same as those returned
    ## by 'sgedgesByGene()' on a SplicingGraphs object of length 1.
    check_way1_vs_way2(
        sgedges(sgi),
        sgedgesByGene(sgi))
    ## After removal of the artificial edges, the edges returned
    ## by 'sgedges( , keep.dup.edges=TRUE)' should be the same as
    ## those returned by 'sgedgesByGene( , keep.dup.edges=TRUE)' or by
    ## 'sgedgesByTranscript()' on a SplicingGraphs object of length 1.
    res1 <- DataFrame(sgedges(sgi, keep.dup.edges=TRUE))
    check_way1_vs_way2(
        res1,
        sgedgesByGene(sgi, keep.dup.edges=TRUE))
    check_way1_vs_way2(
        res1,
        sgedgesByTranscript(sgi))
}
}
