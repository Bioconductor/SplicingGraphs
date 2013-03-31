\name{sgedgesByTranscript-methods}

\alias{sgedgesByTranscript-methods}

\alias{intronsByTranscript,SplicingGraphs-method}

\alias{sgedgesByTranscript}
\alias{sgedgesByTranscript,SplicingGraphs-method}

\alias{sgedgesByGene}
\alias{sgedgesByGene,SplicingGraphs-method}


\title{
  Extract all the introns or splicing graph edges from a SplicingGraphs object
}

\description{
  \code{intronsByTranscript} extracts all the introns of all
  the genes from a \link{SplicingGraphs} object.

  \code{sgedgesByTranscript} extracts all the splicing graph edges of all
  the genes from a \link{SplicingGraphs} object.

  Both functions above return a named \link[GenomicRanges]{GRangesList} object
  where the items (i.e. introns for the former, splicing graph edges for
  the latter) are grouped by transcript. The names on the object are the
  gene ids.

  Alternatively \code{sgedgesByGene} extracts the splicing graph edges like
  \code{sgedgesByTranscript} does but groups them by gene instead of by
  transcript. In addition, within each group, edges with the same
  \emph{global edge id} are merged into a single element.
}

\usage{
\S4method{intronsByTranscript}{SplicingGraphs}(x)

sgedgesByTranscript(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)

sgedgesByGene(x, with.exon.mcols=FALSE, with.hits.mcols=FALSE)
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

    \item \code{\link{countReads}} for assigning reads to a
          \link{SplicingGraphs} object and counting them.
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
## 2. intronsByTranscript()
## ---------------------------------------------------------------------
## The introns of all the genes can be extracted with
## intronsByTranscript(). The result is a *named* GRangesList object
## containing the introns grouped by transcript. The names on the
## object are the gene ids:
in_by_tx <- intronsByTranscript(sg)
in_by_tx

## A sanity check:
ex_by_tx <- unlist(sg)
stopifnot(identical(elementLengths(in_by_tx) + 1L,
                    elementLengths(ex_by_tx)))

## ---------------------------------------------------------------------
## 3. sgedgesByTranscript()
## ---------------------------------------------------------------------
## 'sgedgesByTranscript(sg)' combines 'ex_by_tx' and 'in_by_tx' in
## a single GRangesList object:
edges_by_tx <- sgedgesByTranscript(sg)

## A sanity check:
stopifnot(identical(elementLengths(edges_by_tx),
                    elementLengths(ex_by_tx) + elementLengths(in_by_tx)))

## Select elements that correspond to geneB: 
edges_by_tx[names(edges_by_tx) \%in\% "geneB"]
## Note that edge with global edge id "geneB:3,4" is an intron that
## belongs to transcripts B1 and B2.

## Keep the "exon metadata columns".
sgedgesByTranscript(sg, with.exon.mcols=TRUE)
## Note that those cols are set to NA for edges of type intron.

## ---------------------------------------------------------------------
## 4. sgedgesByGene()
## ---------------------------------------------------------------------
## 'sgedgesByGene' groups the edges by gene instead of by transcript:
edges_by_gene <- sgedgesByGene(sg)
edges_by_gene[["geneB"]]
## Note that the 2 inner elements in 'edges_by_tx' with global edge
## id "geneB:3,4" are now merged into a single inner element.

sgedgesByGene(sg, with.exon.mcols=TRUE)
}
