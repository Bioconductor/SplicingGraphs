\name{rsgedgesByGene-methods}

\alias{rsgedgesByGene-methods}

\alias{uninformativeSSids}
\alias{uninformativeSSids,ANY-method}
\alias{uninformativeSSids,DataFrame-method}

\alias{rsgedgesByTranscript}
\alias{rsgedgesByTranscript,SplicingGraphs-method}

\alias{rsgedgesByGene}
\alias{rsgedgesByGene,SplicingGraphs-method}

\alias{rsgedges}
\alias{sgedges2}

\alias{rsgraph}
\alias{sgraph2}


\title{
  Extract the reduced edges and their ranges from a SplicingGraphs object
}

\description{
  \code{rsgedgesByGene} and \code{rsgedgesByTranscript} are analog to
  \code{\link{sgedgesByGene}} and \code{\link{sgedgesByTranscript}},
  but operate on the \emph{reduced} splicing graphs, that is, the
  graphs in \link{SplicingGraphs} object \code{x} are reduced before
  the edges and their ranges are extracted. The reduced graphs are
  obtained by removing the uninformative nodes from it. See Details
  section below.

  \code{rsgedges} extracts the edges of the reduced splicing graph of
  a given gene from a \link{SplicingGraphs} object.

  \code{rsgraph} extracts the reduced splicing graph for a given gene
  from a \link{SplicingGraphs} object, and returns it as a plottable
  graph-like object.
}

\usage{
rsgedgesByGene(x, with.hits.mcols=FALSE, keep.dup.edges=FALSE)

rsgedgesByTranscript(x, with.hits.mcols=FALSE)

rsgedges(x)

rsgraph(x, tx_id.as.edge.label=FALSE, as.igraph=FALSE)

## Related utility:
uninformativeSSids(x)
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object. Must be of length 1 for \code{rsgedges},
    \code{rsgraph}, and \code{uninformativeSSids}.
  }
  \item{with.hits.mcols}{
    Whether or not to include the \emph{hits metadata columns} in the
    returned object. See \code{?\link{countReads}} for more information.
  }
  \item{keep.dup.edges}{
    Not supported yet.
  }
  \item{tx_id.as.edge.label}{
    Whether or not to use the transcript ids as edge labels.
  }
  \item{as.igraph}{
    TODO
  }
}

\details{
  TODO: Explain graph reduction.
}

\value{
  For \code{rsgedgesByGene}: A \link[GenomicRanges]{GRangesList} object
  named with the gene ids and where the reduced splicing graph edges are
  grouped by gene.

  TODO: Explain values returned by the other function.
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
example(SplicingGraphs)
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

## ---------------------------------------------------------------------
## 2. rsgedgesByGene()
## ---------------------------------------------------------------------
edges_by_gene <- rsgedgesByGene(sg)
edges_by_gene
## 'edges_by_gene' has the length and names of 'sg', that is, the names
## on it are the gene ids and are guaranteed to be unique.

## Extract the reduced edges and their ranges for a given gene:
edges_by_gene[["geneA"]]
## Note that edge with global reduced edge id "geneA:1,2,4,5" is a mixed
## edge obtained by combining together edges "geneA:1,2" (exon),
## "geneA:2,4" (intron), and "geneA:4,5" (exon), during the graph
## reduction.

stopifnot(identical(edges_by_gene["geneB"], rsgedgesByGene(sg["geneB"])))

## ---------------------------------------------------------------------
## 3. sgedgesByTranscript()
## ---------------------------------------------------------------------
#edges_by_tx <- rsgedgesByTranscript(sg)  # not ready yet!
#edges_by_tx

## ---------------------------------------------------------------------
## 4. rsgedges(), rsgraph(), uninformativeSSids()
## ---------------------------------------------------------------------
plot(sgraph(sg["geneB"]))
uninformativeSSids(sg["geneB"])

plot(rsgraph(sg["geneB"]))
rsgedges(sg["geneB"])

## ---------------------------------------------------------------------
## 5. Sanity checks
## ---------------------------------------------------------------------
## TODO: Do the same kind of sanity checks that are done for sgedges()
## vs sgedgesByGene() vs sgedgesByTranscript() (in man page for sgedges).
}
