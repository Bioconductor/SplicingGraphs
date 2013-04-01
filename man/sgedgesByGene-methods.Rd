\name{sgedgesByGene-methods}

\alias{sgedgesByGene-methods}

\alias{intronsByTranscript,SplicingGraphs-method}

\alias{sgedgesByTranscript}
\alias{sgedgesByTranscript,SplicingGraphs-method}

\alias{sgedgesByGene}
\alias{sgedgesByGene,SplicingGraphs-method}


\title{
  Extract all the edges and ranges from a SplicingGraphs object
}

\description{
  \code{sgedgesByGene} and \code{sgedgesByTranscript} both extract the
  edges and ranges of all the genes from a \link{SplicingGraphs} object.
  They return them in a \link[GenomicRanges]{GRangesList} object named
  with the gene ids, and where the items are grouped by gene (for
  \code{sgedgesByGene}) or by transcript (for \code{sgedgesByTranscript}).

  Alternatively, \code{intronsByTranscript} extracts the intronic edges
  and ranges of all the genes from a \link{SplicingGraphs} object.
  It returns them in a \link[GenomicRanges]{GRangesList} object named
  with the gene ids, and where the items are grouped by transcript.
}

\usage{
sgedgesByGene(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE,
                 keep.dup.edges=FALSE)

sgedgesByTranscript(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)

\S4method{intronsByTranscript}{SplicingGraphs}(x)
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object.
  }
  \item{with.exon.mcols}{
    Whether to include or not the \emph{exon metadata columns} in the
    returned object. Those columns are named: \code{exon_id},
    \code{exon_name}, \code{exon_rank}, \code{start_SSid}, and
    \code{end_SSid}. They are set to \code{NA} for edges of type intron.
  }
  \item{with.hits.mcols}{
    TODO
  }
  \item{keep.dup.edges}{
    If \code{FALSE} (the default), then within each group of the returned
    object, edges with the same \emph{global edge id} are merged into
    a single element. Use \code{keep.dup.edges=TRUE} if this merging should
    not be performed.
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

    \item \code{\link{txpath}} for extracting the transcript paths of a
          splicing graph.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{bubbles}} for computing the bubbles of a splicing graph.

    \item \code{\link{countReads}} for assigning reads to the edges of a
          \link{SplicingGraphs} object and summarizing them.
  }
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
## 2. sgedgesByGene()
## ---------------------------------------------------------------------
edges_by_gene <- sgedgesByGene(sg)
edges_by_gene
## 'edges_by_gene' has the length and names of 'sg', that is, the names
## on it are the gene ids and are guaranteed to be unique.

## Extract the edges and ranges for a given gene:
edges_by_gene[["geneB"]]
## Note that edge with global edge id "geneB:3,4" is an intron that
## belongs to transcripts B1 and B2.

edges_by_gene0 <- sgedgesByGene(sg, keep.dup.edges=TRUE)
edges_by_gene0[["geneB"]]
## Note that edge "geneB:3,4" now shows up twice, once for transcript
## B1, and once for transcript B2.

## Keep the "exon metadata columns":
sgedgesByGene(sg, with.exon.mcols=TRUE)
## Note that those cols are set to NA for intronic edges.

## ---------------------------------------------------------------------
## 3. sgedgesByTranscript()
## ---------------------------------------------------------------------
edges_by_tx <- sgedgesByTranscript(sg)
edges_by_tx

## 'edges_by_tx' is typically longer than 'sg'.
## IMPORTANT NOTE: One caveat here is that the names on 'edges_by_tx'
## are the gene ids, not the transcript ids, and thus are typically NOT
## unique!

## Select elements of a given gene:
edges_by_tx["geneB"]  # not a good idea
edges_by_tx[names(edges_by_tx) \%in\% "geneB"]  # much better :-)
## Note that edge with global edge id "geneB:3,4" is an intron that
## belongs to transcripts B1 and B2.

## Keep the "exon metadata columns":
sgedgesByTranscript(sg, with.exon.mcols=TRUE)
## Note that those cols are set to NA for intronic edges.

## ---------------------------------------------------------------------
## 4. intronsByTranscript()
## ---------------------------------------------------------------------
in_by_tx <- intronsByTranscript(sg)
in_by_tx

## 'in_by_tx' has the length and names of 'edges_by_tx'. The same
## recommendation applies for selecting elements of a given set of
## genes:
in_by_tx[c("geneB", "geneD")]  # not a good idea
in_by_tx[names(in_by_tx) \%in\% c("geneB", "geneD")]  # much better :-)

## ---------------------------------------------------------------------
## 5. Comparing the outputs of unlist(), intronsByTranscript(), and
##    sgedgesByTranscript()
## ---------------------------------------------------------------------
ex_by_tx <- unlist(sg)
in_by_tx <- intronsByTranscript(sg)
edges_by_tx <- sgedgesByTranscript(sg)

## A sanity check:
stopifnot(identical(elementLengths(in_by_tx) + 1L,
                    elementLengths(ex_by_tx)))

## 'edges_by_tx' combines 'ex_by_tx' and 'in_by_tx' in a single
## GRangesList object. Sanity check:
stopifnot(identical(elementLengths(edges_by_tx),
                    elementLengths(ex_by_tx) + elementLengths(in_by_tx)))
}
