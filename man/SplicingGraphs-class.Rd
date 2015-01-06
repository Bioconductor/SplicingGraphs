\name{SplicingGraphs-class}

\alias{class:GeneModel}
\alias{GeneModel-class}

\alias{class:SplicingGraphs}
\alias{SplicingGraphs-class}

\alias{seqnames,GeneModel-method}
\alias{strand,GeneModel-method}
\alias{seqinfo,GeneModel-method}

\alias{length,SplicingGraphs-method}
\alias{names,SplicingGraphs-method}
\alias{[,SplicingGraphs,ANY,ANY-method}
\alias{[[,SplicingGraphs,ANY,ANY-method}
\alias{elementLengths,SplicingGraphs-method}
\alias{unlist,SplicingGraphs-method}
\alias{seqnames,SplicingGraphs-method}
\alias{strand,SplicingGraphs-method}
\alias{seqinfo,SplicingGraphs-method}

\alias{show,SplicingGraphs-method}

\alias{SplicingGraphs}
\alias{SplicingGraphs,GRangesList-method}
\alias{SplicingGraphs,TxDb-method}


\title{
  SplicingGraphs objects
}

\description{
  The SplicingGraphs class is a container for storing splicing graphs
  together with the gene model that they are based on.
}

\usage{
SplicingGraphs(x, grouping=NULL, min.ntx=2, max.ntx=NA, check.introns=TRUE)

## SplicingGraphs basic API:

\S4method{length}{SplicingGraphs}(x)

\S4method{names}{SplicingGraphs}(x)

\S4method{seqnames}{SplicingGraphs}(x)

\S4method{strand}{SplicingGraphs}(x)

\S4method{[}{SplicingGraphs,ANY,ANY}(x, i, j, ... , drop=TRUE)

\S4method{[[}{SplicingGraphs,ANY,ANY}(x, i, j, ...)

\S4method{elementLengths}{SplicingGraphs}(x)

\S4method{unlist}{SplicingGraphs}(x, recursive=TRUE, use.names=TRUE)

\S4method{seqinfo}{SplicingGraphs}(x)
}

\arguments{
  \item{x}{
    For \code{SplicingGraphs}: A \link[GenomicRanges]{GRangesList} object
    containing the exons of one or more genes grouped by transcript.
    Alternatively, \code{x} can be a \link[GenomicFeatures]{TxDb}
    object. See Details section below.

    For the methods in the SplicingGraphs basic API: A SplicingGraphs object.
  }
  \item{grouping}{
    An optional object that represents the grouping by gene of the top-level
    elements (i.e. the transcripts) in \code{x}. See Details section
    below.
  }
  \item{min.ntx, max.ntx}{
    Single integers (or \code{NA} for \code{max.ntx}) specifying the minimum
    and maximum number of transcripts a gene must have to be considered for
    inclusion in the object returned by \code{SplicingGraphs}. A value of
    \code{NA} for \code{max.ntx} means no maximum.
  }
  \item{check.introns}{
    If \code{TRUE}, \code{SplicingGraphs} checks that, within each transcript,
    exons are ordered from 5' to 3' with gaps of at least 1 nucleotide between
    them.
  }
  \item{i, j, ..., drop}{
    A SplicingGraphs object is a list-like object and therefore it can be
    subsetted like a list. When subsetting with \code{[}, the result is
    another SplicingGraphs object containing only the selected genes.
    When subsetting with \code{[[}, the result is an \emph{unnamed}
    \link[GenomicRanges]{GRangesList} object containing the exons grouped
    by transcript.
    Like for list, subsetting only accepts 1 argument (\code{i}). The
    \code{drop} argument is ignored and trying to pass any additional
    argument (to \code{j} or in \code{...}) will raise an error.
  }
  \item{recursive, use.names}{
    A SplicingGraphs object is a list-like object and therefore it can be
    unlisted with \code{unlist}. The result is a
    \link[GenomicRanges]{GRangesList} object containing the exons
    grouped by transcript. By default this object has names on it, and
    the names are the gene ids. Note that because each element in this
    object represents a transcript (and not a gene), the names are not
    unique. If \code{use.names=FALSE} is used, the result has no names on it.
    The \code{recursive} agument is ignored.
  }
}

\details{
  The Splicing graph theory only applies to genes that have all the exons
  of all their transcripts on the same chromosome and strand. In particular,
  in its current form, the splicing graph theory cannot describe
  trans-splicing events. The \code{SplicingGraphs} constructor will reject
  genes that do not satisfy this.

  The first argument of the \code{SplicingGraphs} constructor, \code{x},
  can be either a \link[GenomicRanges]{GRangesList} object or a
  \link[GenomicFeatures]{TxDb} object.

  When \code{x} is a \link[GenomicRanges]{GRangesList} object, it must
  contain the exons of one or more genes grouped by transcripts. More
  precisely, each top-level element in \code{x} must contain the genomic
  ranges of the exons for a particular transcript. Typically \code{x}
  will be obtained from a \link[GenomicFeatures]{TxDb} object
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
    \item A data.frame or \link[S4Vectors]{DataFrame} with 2 character vector
          columns: a gene_id column (factor, character vector, or integer
          vector), and a tx_id (or alternatively tx_name) column. In that
          case, \code{x} must be named and all the values in the tx_id
          (or tx_name) column must be in \code{names(x)}.
  }
}

\value{
  For \code{SplicingGraphs}: a SplicingGraphs object with 1 element per gene.

  For \code{length}: the number of genes in \code{x}, which is also the number
  of splicing graphs in \code{x}.

  For \code{names}: the gene ids. Note that the names on a SplicingGraphs
  object are always unique and cannot be modified.

  For \code{seqnames}: a named factor of the length of \code{x} containing
  the name of the chromosome for each gene. 

  For \code{strand}: a named factor of the length of \code{x} containing
  the strand for each gene.

  For \code{elementLengths}: the number of transcripts per gene.

  For \code{seqinfo}: the seqinfo of the \link[GenomicRanges]{GRangesList}
  or \link[GenomicFeatures]{TxDb} object that was used to construct
  the SplicingGraphs object.
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
  This man page is part of the \pkg{SplicingGraphs} package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.

  Other topics related to this man page and documented in other packages:
  \itemize{
    \item The \code{\link[GenomicFeatures]{exonsBy}} and
          \code{\link[GenomicFeatures]{transcriptsBy}} functions, and
          the \link[GenomicFeatures]{TxDb} class, defined in the
          \pkg{GenomicFeatures} package.

    \item The \link[GenomicRanges]{GRangesList} class defined in the
          \pkg{GenomicRanges} package.

    \item The \link[IRanges]{IntegerList} and \link[IRanges]{CharacterList}
          classes defined in the \pkg{IRanges} package.

    \item The \link[S4Vectors]{DataFrame} class defined in the
          \pkg{S4Vectors} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. Load a toy gene model as a TxDb object
## ---------------------------------------------------------------------

library(GenomicFeatures)
suppressWarnings(
  toy_genes_txdb <- makeTxDbFromGFF(toy_genes_gff())
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

## Alternatively 'sg' can be constructed directly from the TxDb
## object:
sg2 <- SplicingGraphs(toy_genes_txdb)  # same as 'sg'
sg2

## Note that because SplicingGraphs objects have a slot that is an
## environment (for caching the bubbles), they cannot be compared with
## 'identical()' (will always return FALSE). 'all.equal()' should be
## used instead:
stopifnot(isTRUE(all.equal(sg2, sg)))

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids:
length(sg)
names(sg)

## ---------------------------------------------------------------------
## 3. Basic manipulation of a SplicingGraphs object
## ---------------------------------------------------------------------

## Basic accessors:
seqnames(sg)
strand(sg)
seqinfo(sg)

## Number of transcripts per gene:
elementLengths(sg)

## The transcripts of a given gene can be extracted with [[. The result
## is an *unnamed* GRangesList object containing the exons grouped by
## transcript:
sg[["geneD"]]

## See '?plotTranscripts' for how to plot those transcripts.

## The transcripts of all the genes can be extracted with unlist(). The
## result is a *named* GRangesList object containing the exons grouped
## by transcript. The names on the object are the gene ids:
ex_by_tx <- unlist(sg)
ex_by_tx
}
