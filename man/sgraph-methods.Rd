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
sgraph(x, gene_id=NA, keep.dup.edges=FALSE,
       tx_id.as.edge.label=FALSE, as.igraph=FALSE)
sgraph2(x, gene_id=NA,
        tx_id.as.edge.label=FALSE, as.igraph=FALSE)

## PLotting:

\S4method{plot}{SplicingGraphs,ANY}(x, y, gene_id=NA)
slideshow(x)
}

\arguments{
  \item{x}{
    TODO
  }
  \item{gene_id}{
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

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{bubbles}} for extracting the bubbles of a splicing graph.
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

sgA <- sgraph(sg, gene_id="geneA", tx_id.as.edge.label=TRUE)

if (interactive()) {
  ## Edges are labeled with the transcript ids (or names), in blue.
  ## The orange arrows are edges corresponding to exons:
  plot(sgA)

  ## Note that plot() works directly on the SplicingGraphs object 'sg':
  plot(sg, "geneB")

  ## Slideshow of the graphs:
  slideshow(sg)
}
}
