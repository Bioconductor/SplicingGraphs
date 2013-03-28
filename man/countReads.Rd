\name{countReads}

\alias{assignReads}
\alias{countReads}

\title{
  Assign reads to splicing graphs and count them
}

\description{
  \code{assignReads} assigns reads to the exons and introns of a
  \link{SplicingGraphs} object.
  \code{countReads} count the assigned reads.
}

\usage{
assignReads(sg, reads, sample.name=NA)
countReads(sg)
}

\arguments{
  \item{sg}{
    A \link{SplicingGraphs} object.
  }
  \item{reads}{
    A \link[GenomicRanges]{GRangesList} object containing the reads
    to assign to the exons and introns in \code{sg}.
  }
  \item{sample.name}{
    A single string containing the name of the sample where the reads
    are coming from.
  }
}

\details{
  TODO
}

\value{
  For \code{assignReads}: the supplied \link{SplicingGraphs} object with
  the reads assigned to it.

  For \code{countReads}: a \link[IRanges]{DataFrame} object with one row
  per unique splicing graph edge and one column of counts per sample.
  Two additional columns (\code{"sgedge_id"} and \code{"ex_or_in"}) contain
  the splicing graph edge ids and the type of edge (exon or intron),
  respectively.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \link{SplicingGraphs} class.

    \item \code{\link{sgedgesByTranscript}} for extracting all the introns
          or splicing graph edges from a \link{SplicingGraphs} object.

    \item \code{\link{txpath}} for extracting the transcript paths of a
          splicing graph.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{bubbles}} for computing the bubbles of a splicing graph.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. Make SplicingGraphs object 'sg' from toy gene model (see
##    '?SplicingGraphs')
## ---------------------------------------------------------------------

example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

## ---------------------------------------------------------------------
## 2. Load toy reads
## ---------------------------------------------------------------------

## Load toy reads (single-end) from a BAM file. We filter out secondary
## alignments, reads not passing quality controls, and PCR or optical
## duplicates (see ?scanBamFlag for more information):
flag0 <- scanBamFlag(isNotPrimaryRead=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isDuplicate=FALSE)
param0 <- ScanBamParam(flag=flag0)
gal <- readGappedAlignments(toy_reads_bam(), use.names=TRUE, param=param0)
gal

## Put the reads in a GRangesList object. Using 'order.as.in.query=TRUE'
## is important to ensure that the ranges corresponding to each read are
## in the correct order.
grl <- grglist(gal, order.as.in.query=TRUE)

## ---------------------------------------------------------------------
## 3. Assign the reads to the exons and introns in 'sg'
## ---------------------------------------------------------------------
## The same read can be assigned to more than 1 exon or intron (e.g. a
## junction read with 1 gap can be assigned to 2 exons and 1 intron).
sg <- assignReads(sg, grl)

sg[["geneA"]]
sgedges(sg["geneA"])

## ---------------------------------------------------------------------
## 4. Count the number of reads per splicing graph edge
## ---------------------------------------------------------------------
countReads(sg)

## ---------------------------------------------------------------------
## 5. Count nb of unambiguous compatible hits per transcript
## ---------------------------------------------------------------------

#ov2 <- ov1[queryHits(ov1) \%in\% which(countQueryHits(ov1) == 1L)]
#UATXHcount <- countSubjectHits(ov2)
#mcols(sg@genes@unlistData)$UATXHcount <- UATXHcount

#plot(sg["geneA"])
}
