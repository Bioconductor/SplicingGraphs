\name{SplicingGraphs-class}

\alias{class:.SplicingGraphGenes}
\alias{.SplicingGraphGenes-class}

\alias{class:SplicingGraphs}
\alias{SplicingGraphs-class}

\alias{SplicingGraphs}
\alias{SplicingGraphs,GRangesList-method}
\alias{SplicingGraphs,TranscriptDb-method}

\alias{length,SplicingGraphs-method}
\alias{names,SplicingGraphs-method}
\alias{[,SplicingGraphs,ANY,ANY-method}
\alias{[[,SplicingGraphs,ANY,ANY-method}
\alias{elementLengths,SplicingGraphs-method}
\alias{unlist,SplicingGraphs-method}
\alias{intronsByTranscript,SplicingGraphs-method}

\alias{show,SplicingGraphs-method}

\alias{plotTranscripts}
\alias{plotTranscripts,GRangesList-method}
\alias{plotTranscripts,TranscriptDb-method}
\alias{plotTranscripts,SplicingGraphs-method}


\title{
  SplicingGraphs objects
}

\description{
  The SplicingGraphs class is a container for splicing graphs created from
  transcript annotations.
}

\usage{
SplicingGraphs(x, grouping=NULL, min.ntx=2, max.ntx=NA, check.introns=TRUE)

## SplicingGraphs basic API:

\S4method{length}{SplicingGraphs}(x)

\S4method{names}{SplicingGraphs}(x)

\S4method{[}{SplicingGraphs,ANY,ANY}(x, i, j, ... , drop=TRUE)

\S4method{[[}{SplicingGraphs,ANY,ANY}(x, i, j, ...)

\S4method{elementLengths}{SplicingGraphs}(x)

\S4method{unlist}{SplicingGraphs}(x, recursive=TRUE, use.names=TRUE)

\S4method{intronsByTranscript}{SplicingGraphs}(x)

plotTranscripts(x)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRangesList} object containing the exons of one
    or more genes grouped by transcripts. Alternatively, \code{x} can be a
    \link[GenomicFeatures]{TranscriptDb} object. See Details section below.
  }
  \item{grouping}{
    An optional object that represents the grouping by gene of the top-level
    elements (i.e. the transcripts) in \code{x}. See Details section
    below.
  }
  \item{min.ntx}{
    TODO
  }
  \item{max.ntx}{
    TODO
  }
  \item{check.introns}{
    TODO
  }
  \item{i, j, ..., drop}{
    TODO
  }
  \item{recursive, use.names}{
    TODO
  }
}

\details{
  The first argument of the \code{SplicingGraphs} constructor, \code{x},
  can be either a \link[GenomicRanges]{GRangesList} object or a
  \link[GenomicFeatures]{TranscriptDb} object.

  When \code{x} is a \link[GenomicRanges]{GRangesList} object, it must
  contain the exons of one or more genes grouped by transcripts. More
  precisely, each top-level element in \code{x} must contain the genomic
  ranges of the exons for a particular transcript. Typically \code{x}
  will be obtained from a \link[GenomicFeatures]{TranscriptDb} object
  \code{txdb} with \code{\link[GenomicFeatures]{exonsBy}(txdb, by="tx",
  use.names=TRUE)}.

  \code{grouping} is an optional argument that is only supported
  when \code{x} is a \link[GenomicRanges]{GRangesList} object.
  It represents the grouping by gene of the top-level elements (i.e. the
  transcripts) in \link[GenomicRanges]{GRangesList} object \code{x}.
  It can be either:
  \itemize{
    \item Missing (i.e. \code{NULL}). In that case, all the transcripts in
          \code{x} are considered to belong to the same gene and the
          SplicingGraphs object returned by \code{SplicingGraphs} will be
          unnamed.
    \item A list of integer or character vectors, or an
          \link[IRanges]{IntegerList}, or a \link[IRanges]{CharacterList}
          object, of length the number of genes to process, and where
          \code{grouping[[i]]} is a vector of valid subscripts in
          \code{x} pointing to all the transcripts of the i-th gene.
    \item A factor, character vector, or integer vector, of the same length
          as \code{x} and 1 level per gene.
    \item A named \link[GenomicRanges]{GRangesList} object containing
          transcripts grouped by genes i.e. each top-level element in
          \code{grouping} contains the genomic ranges of the transcripts
          for a particular gene. In that case, the grouping is inferred from
          the tx_id (or alternatively tx_name) metadata column of
          \code{unlist(grouping)} and all the values in that column must
          be in \code{names(x)}.
          If \code{x} was obtained with
          \code{\link[GenomicFeatures]{exonsBy}(txdb, by="tx", use.names=TRUE)},
          then the \link[GenomicRanges]{GRangesList} object used for grouping
          would typically be obtained with
          \code{\link[GenomicFeatures]{transcriptsBy}(txdb, by="gene")}.
    \item A data.frame or \link[IRanges]{DataFrame} with 2 character vector
          columns: a gene_id column (factor, character vector, or integer
          vector), and a tx_id (or alternatively tx_name) column. In that
          case, \code{x} must be named and all the values in the tx_id
          (or tx_name) column must be in \code{names(x)}.
  }
}

\value{
  TODO
}

\author{
  H. Pages
}

\references{
  Heber, S., Alekseyev, M., Sze, S., Tang, H., and Pevzner, P. A. 
  \emph{Splicing graphs and EST assembly problem}
  Bioinformatics
  Date: Jul 2002
  Vol: 18
  Pages: S181-S188

  Sammeth, M. (2009) 
  \emph{Complete Alternative Splicing Events Are Bubbles in Splicing Graphs}
  J. Comput. Biol.
  Date: Aug 2009
  Vol: 16
  Pages: 1117-1140
}

\seealso{
  \itemize{
    \item \code{\link[GenomicFeatures]{exonsBy}},
          \code{\link[GenomicFeatures]{transcriptsBy}}, and the
          \link[GenomicFeatures]{TranscriptDb} class in the GenomicFeatures
          package.

    \item The \link[GenomicRanges]{GRangesList} class in the GenomicRanges
          package.

    \item The \link[IRanges]{IntegerList}, \link[IRanges]{CharacterList},
          and \link[IRanges]{DataFrame} classes in the IRanges package.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a plottable
          graph-like object.

    \item \code{\link{bubbles}} for extracting the bubbles of a splicing graph.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. Load a toy gene model as a TranscriptDb object
## ---------------------------------------------------------------------

library(GenomicFeatures)
suppressWarnings(
  toy_genes_txdb <- makeTranscriptDbFromGFF(toy_genes_gff())
)

## ---------------------------------------------------------------------
## 2. Compute all the splicing graphs (1 graph per gene) and return them
##    in a SplicingGraphs object
## ---------------------------------------------------------------------

## Extract the exons grouped by transcript:
ex_by_tx <- exonsBy(toy_genes_txdb, by="tx", use.names=TRUE)

## Extract the transcripts grouped by gene:
tx_by_gn <- transcriptsBy(toy_genes_txdb, by="gene")

sg <- SplicingGraphs(ex_by_tx, tx_by_gn)
sg

## Alternatively 'sg' can be constructed directly from the TranscriptDb
## object:
sg2 <- SplicingGraphs(toy_genes_txdb)  # same as 'sg'
sg2

## Note that because SplicingGraphs objects have a slot that is an
## environment (for caching the bubbles), they cannot be compared with
## identical() (will always return FALSE). A workaround is to unlist
## them before comparison (see below for what unlist() does exactly on
## a SplicingGraphs object):
stopifnot(identical(unlist(sg2), unlist(sg)))

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids:
length(sg)
names(sg)

## The transcripts of a given gene can be extracted with [[. The result
## is an *unnamed* GRangesList object containing the exons grouped by
## transcript:
sg[["geneD"]]
plotTranscripts(sg[["geneD"]])

## The transcripts of all the genes can be extracted with unlist(). The
## result is a *named* GRangesList object containing the exons grouped
## by transcript. The names of the object are the gene ids:
ex_by_tx <- unlist(sg)
ex_by_tx
plotTranscripts(ex_by_tx)

## Similarly, the introns can be extracted with intronsByTranscript().
## The result is a *named* GRangesList object containing the introns
## grouped by transcript. The names of the object are the gene ids:
in_by_tx <- intronsByTranscript(sg)
in_by_tx

## A sanity check:
stopifnot(identical(elementLengths(in_by_tx) + 1L, elementLengths(ex_by_tx)))
}
