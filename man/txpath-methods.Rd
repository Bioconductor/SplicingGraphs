\name{txpath-methods}

\alias{txpath-methods}

\alias{txpath}
\alias{txpath,GRangesList-method}
\alias{txpath,SplicingGraphs-method}

\alias{txweight}
\alias{txweight,SplicingGraphs-method}
\alias{txweight<-}
\alias{txweight<-,SplicingGraphs-method}


\title{
  Extract the transcript paths of a splicing graph
}

\description{
  \code{txpath} extracts the transcript paths of the splicing graph
  of a given gene from a \link{SplicingGraphs} object.
}

\usage{
txpath(x, as.matrix=FALSE)

## Related utility:
txweight(x)
txweight(x) <- value
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object of length 1
    or a \link[GenomicRanges]{GRangesList} object for \code{txpath}.

    A \link{SplicingGraphs} object for \code{txweight}.
  }
  \item{as.matrix}{
    TODO
  }
  \item{value}{
    A numeric vector containing the weights to assign to each
    transcript in \code{x}.
  }
}

\details{
  TODO
}

\value{
  A named list-like object with one list element per transcript in the gene.
  Each list element is an integer vector that describes the \emph{path}
  of the transcript i.e. the \emph{Splicing Site ids} that it goes thru.
}

\author{
  H. Pages
}

\seealso{
  This man page is part of the SplicingGraphs package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.

  Other topics related to this man page and documented in other packages:
  \itemize{
    \item The \link[GenomicRanges]{GRangesList},
          \link[GenomicRanges]{GappedAlignments}, and
          \link[GenomicRanges]{GappedAlignmentPairs} classes
          in the GenomicRanges package.

    \item \link[GenomicRanges]{findOverlaps-methods} and
          \link[GenomicRanges]{encodeOverlaps-methods} in the
          GenomicRanges package.

    \item The \code{\link[Rsamtools]{ScanBamParam}} function in the
          Rsamtools package.
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
## 2. txpath()
## ---------------------------------------------------------------------
## Note that the list elements in the returned IntegerList object
## always consist of an even number of Splicing Site ids in ascending
## order.
txpath(sg["geneB"])
txpath(sg["geneD"])
strand(sg)

txpath(sg["geneD"], as.matrix=TRUE)  # splicing matrix

## ---------------------------------------------------------------------
## 3. txweight()
## ---------------------------------------------------------------------
txweight(sg)
plot(sg["geneD"])

txweight(sg) <- 5
txweight(sg)
plot(sg["geneD"])  # FIXME: Edges not rendered with correct width!
plot(sgraph(sg["geneD"], as.igraph=TRUE))  # need to use this for now

txweight(sg)[8:11] <- 5 * (4:1)
txweight(sg)
plot(sgraph(sg["geneD"], tx_id.as.edge.label=TRUE, as.igraph=TRUE))

## ---------------------------------------------------------------------
## 4. An advanced example
## ---------------------------------------------------------------------
## [TODO: Counting "unambiguous compatible hits" per transcript should be
##  supported by countReads(). Simplify the code below when countReads()
##  supports this.]
## Here we show how to find "unambiguous compatible hits" between a set
## of RNA-seq reads and a set of transcripts, that is, hits that are
## compatible with the splicing of exactly 1 transcript. Then we set the
## transcript weights based on the number of unambiguous compatible hits
## they received and we finally plot some splicing graphs that show
## the weighted transcripts.
## Note that the code we use to find the unambiguous compatible hits
## uses findOverlaps() and encodeOverlaps() defined in the IRanges and
## GenomicRanges packages. It only requires that the transcripts are
## represented as a GRangesList object and the reads as a
## GappedAlignments (single-end) or GappedAlignmentPairs (paired-end)
## object, and therefore is not specific to SplicingGraphs.

## First we load toy reads (single-end) from a BAM file. We filter out
## secondary alignments, reads not passing quality controls, and PCR or
## optical duplicates (see ?scanBamFlag in the Rsamtools package for
## more information):
flag0 <- scanBamFlag(isNotPrimaryRead=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isDuplicate=FALSE)
param0 <- ScanBamParam(flag=flag0)
gal <- readGappedAlignments(toy_reads_bam(), use.names=TRUE, param=param0)
gal

## Put the reads in a GRangesList object:
grl <- grglist(gal, order.as.in.query=TRUE)

## Put the transcripts in a GRangesList object (made of exons grouped
## by transcript):
ex_by_tx <- unlist(sg)

## Most of the times the RNA-seq protocol is unstranded so the strand
## reported in the BAM file (and propagated to 'grl') for each alignment
## is meaningless. Thus we should call findOverlaps() with
## 'ignore.strand=TRUE':
ov0 <- findOverlaps(grl, ex_by_tx, ignore.strand=TRUE)

## Encode the overlaps (this compare the fragmentation of the reads with
## the splicing of the transcripts):
ovenc <- encodeOverlaps(grl, ex_by_tx, hits=ov0,
                        flip.query.if.wrong.strand=TRUE)
is_comp <- isCompatibleWithSplicing(ovenc)

## Keep compatible overlaps only:
ov1 <- ov0[is_comp]

## Only keep overlaps that are compatible with exactly 1 transcript:
ov2 <- ov1[queryHits(ov1) \%in\% which(countQueryHits(ov1) == 1L)]
nhit_per_tx <- countSubjectHits(ov2)
names(nhit_per_tx) <- names(txweight(sg))
nhit_per_tx

txweight(sg) <- 2 * nhit_per_tx
plot(sgraph(sg["geneA"], tx_id.as.edge.label=TRUE, as.igraph=TRUE))
plot(sgraph(sg["geneB"], tx_id.as.edge.label=TRUE, as.igraph=TRUE))
}
