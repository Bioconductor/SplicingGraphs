\name{bubbles-methods}

\alias{bubbles-methods}

\alias{bubbles}
\alias{bubbles,SplicingGraphs-method}
\alias{bubbles,ANY-method}
\alias{bubbles,IntegerList-method}

\alias{ASCODE2DESC}


\title{
  Compute the bubbles of a splicing graph
}

\description{
  \code{bubbles} computes the bubbles of the splicing graph of a given gene
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
  This man page is part of the \pkg{SplicingGraphs} package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.
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
