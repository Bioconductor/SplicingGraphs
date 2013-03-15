\name{bubbles-methods}

\alias{bubbles-methods}

\alias{bubbles}
\alias{bubbles,SplicingGraphs-method}
\alias{bubbles,ANY-method}
\alias{bubbles,IntegerList-method}

\alias{ASCODE2DESC}


\title{
  Extract the bubbles of a splicing graph
}

\description{
  \code{bubbles} extracts the bubbles of the splicing graph of a given gene
  from a \link{SplicingGraphs} object.
}

\usage{
bubbles(x)

ASCODE2DESC
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object of length 1.
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
  }
}

\examples{
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

plot(sgraph(sg["geneA"], tx_id.as.edge.label=TRUE))
bubbles(sg["geneA"])

plot(sgraph(sg["geneB"], tx_id.as.edge.label=TRUE))
bubbles(sg["geneB"])

plot(sgraph(sg["geneD"], tx_id.as.edge.label=TRUE))
bubbles(sg["geneD"])
}
