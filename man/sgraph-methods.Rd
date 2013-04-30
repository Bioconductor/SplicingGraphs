\name{sgraph-methods}

\alias{sgraph-methods}

\alias{sgraph}
\alias{sgraph,ANY-method}
\alias{sgraph,data.frame-method}
\alias{sgraph,DataFrame-method}
\alias{sgraph,igraph-method}

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

## PLotting:

\S4method{plot}{SplicingGraphs,ANY}(x, y, ...)
slideshow(x)
}

\arguments{
  \item{x}{
    TODO
  }
  \item{keep.dup.edges}{
    If \code{FALSE} (the default), then edges with the same \emph{global
    edge id} are merged together. Use \code{keep.dup.edges=TRUE} if this
    merging should not be performed.
  }
  \item{tx_id.as.edge.label}{
    Whether or not to use the transcript ids as edge labels.
  }
  \item{as.igraph}{
    TODO
  }
  \item{y}{
    TODO
  }
  \item{...}{
    Additional arguments. \code{plot} passes these to
    \code{plot,Ragraph,ANY-method} for use in, e.g., adding a
    \code{main} title to the plot.
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
  This man page is part of the \pkg{SplicingGraphs} package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.
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
