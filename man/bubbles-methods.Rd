\name{bubbles-methods}

\alias{bubbles-methods}

\alias{bubbles}
\alias{bubbles,ANY-method}
\alias{bubbles,IntegerList-method}


\title{
  Extract the bubbles of a splicing graph
}

\description{
  \code{bubbles} extracts the bubbles of the splicing graph of a given gene
  from a \link{SplicingGraphs} object.
}

\usage{
bubbles(x, gene_id=NA)
}

\arguments{
  \item{x}{
    TODO
  }
  \item{gene_id}{
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

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.
  }
}

\examples{
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per transcript, and each transcript is
## assigned a name that is the id of the gene it belongs to. All the
## transcripts belonging to the same gene are guaranteed to be
## consecutive elements in 'sg'.
names(sg)

plot(sgraph(sg, gene_id="geneA", tx_id.as.edge.label=TRUE))
bubbles(sg, gene_id="geneA")

plot(sgraph(sg, gene_id="geneB", tx_id.as.edge.label=TRUE))
bubbles(sg, gene_id="geneB")

plot(sgraph(sg, gene_id="geneD", tx_id.as.edge.label=TRUE))
bubbles(sg, gene_id="geneD")
}
