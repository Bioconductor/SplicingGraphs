\name{sgraph-methods}

\alias{sgraph-methods}

\alias{sgraph}
\alias{sgraph,ANY-method}
\alias{sgraph,data.frame-method}
\alias{sgraph,DataFrame-method}
\alias{sgraph,igraph-method}

\alias{sgraph2}

\alias{plot,SplicingGraphs,ANY-method}
\alias{slideshow}


\title{
  Extract a splicing graph as a plottable graph-like object
}

\description{
  Extract the splicing graph for a given gene from a \link{SplicingGraphs}
  object and return it as a plottable graph-like object.
}

\usage{
sgraph(x, keep.dup.edges=FALSE, tx_id.as.edge.label=FALSE, as.igraph=FALSE)
sgraph2(x, tx_id.as.edge.label=FALSE, as.igraph=FALSE)

## PLotting:

\S4method{plot}{SplicingGraphs,ANY}(x, y)
slideshow(x)
}

\arguments{
  \item{x}{
    TODO
  }
  \item{keep.dup.edges}{
    TODO
  }
  \item{tx_id.as.edge.label}{
    TODO
  }
  \item{as.igraph}{
    TODO
  }
  \item{y}{
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

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

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

graphA <- sgraph(sg["geneA"], tx_id.as.edge.label=TRUE)

if (interactive()) {
  ## Edges are labeled with the transcript ids (or names), in blue.
  ## The orange arrows are edges corresponding to exons:
  plot(graphA)

  ## Note that plot() works directly on a SplicingGraphs object of
  ## length 1:
  plot(sg["geneB"])

  ## Slideshow of the graphs:
  slideshow(sg)
}
}
