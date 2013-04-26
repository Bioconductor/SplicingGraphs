\name{sgedgesByGene-methods}

\alias{sgedgesByGene-methods}

\alias{intronsByTranscript,SplicingGraphs-method}

\alias{sgedgesByTranscript}
\alias{sgedgesByTranscript,SplicingGraphs-method}

\alias{sgedgesByGene}
\alias{sgedgesByGene,SplicingGraphs-method}


\title{
  Extract the edges and their ranges from a SplicingGraphs object
}

\description{
  \code{sgedgesByGene} and \code{sgedgesByTranscript} both extract the
  edges and their ranges of all the genes from a \link{SplicingGraphs} object.
  They return them in a \link[GenomicRanges]{GRangesList} object named
  with the gene ids, and where the items are grouped by gene (for
  \code{sgedgesByGene}) or by transcript (for \code{sgedgesByTranscript}).

  Alternatively, \code{intronsByTranscript} extracts the intronic edges
  and their ranges of all the genes from a \link{SplicingGraphs} object.
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
    Whether or not to include the \emph{exon metadata columns} in the
    returned object. Those columns are named: \code{exon_id},
    \code{exon_name}, \code{exon_rank}, \code{start_SSid}, and
    \code{end_SSid}. They are set to \code{NA} for edges of type intron.
  }
  \item{with.hits.mcols}{
    Whether or not to include the \emph{hits metadata columns} in the
    returned object. See \code{?\link{countReads}} for more information.
  }
  \item{keep.dup.edges}{
    If \code{FALSE} (the default), then within each group of the returned
    object, edges with the same \emph{global edge id} are merged into
    a single element. Use \code{keep.dup.edges=TRUE} if this merging should
    not be performed.
  }
}

\value{
  A \link[GenomicRanges]{GRangesList} object named with the gene ids and
  where the items are grouped by gene (for \code{sgedgesByGene}), or by
  transcript (for \code{sgedgesByTranscript} and \code{intronsByTranscript}).
  In the latter case (i.e. grouping by transcript), the names are not unique.

  The items that are being grouped are the splicing graph edges of type
  exon and intron (no artificial edges) for \code{sgedgesByGene} and
  \code{sgedgesByTranscript}, and the introns for \code{intronsByTranscript}.

  When the grouping is by transcript (i.e. for \code{sgedgesByTranscript}
  and \code{intronsByTranscript}, items are ordered by their position
  from 5' to 3'.

  About duplicated edges: A given edge can typically be shared by more than
  1 transcript within the same gene, therefore \code{sgedgesByTranscript}
  typically returns an object where the same \emph{global edge id} shows
  up in more than 1 group. However, the same \emph{global edge id} is never
  shared across genes. By default \code{sgedgesByGene} removes duplicated
  edges, unless \code{keep.dup.edges=TRUE} is used.
}

\author{
  H. Pages
}

\seealso{
  This man page is part of the SplicingGraphs package.
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
## 2. sgedgesByGene()
## ---------------------------------------------------------------------
edges_by_gene <- sgedgesByGene(sg)
edges_by_gene
## 'edges_by_gene' has the length and names of 'sg', that is, the names
## on it are the gene ids and are guaranteed to be unique.

## Extract the edges and their ranges for a given gene:
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
