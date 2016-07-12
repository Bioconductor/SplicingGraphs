\name{assignReads}

\alias{assignReads}
\alias{removeReads}


\title{
  Assign reads to the edges of a SplicingGraphs object
}

\description{
  \code{assignReads} assigns reads to the exonic and intronic edges of a
  \link{SplicingGraphs} object.

  \code{removeReads} removes all the reads assigned to the exonic and
  intronic edges of a \link{SplicingGraphs} object.
}

\usage{
assignReads(sg, reads, sample.name=NA)

removeReads(sg)
}

\arguments{
  \item{sg}{
    A \link{SplicingGraphs} object.
  }
  \item{reads}{
    A \link[GenomicAlignments]{GAlignments},
    \link[GenomicAlignments]{GAlignmentPairs}, or
    \link[GenomicRanges]{GRangesList} object, containing the
    reads to assign to the exons and introns in \code{sg}.
    It must have unique names on it, typically the QNAME ("query name")
    field coming from the BAM file. More on this in the 'About the read
    names' section below.
  }
  \item{sample.name}{
    A single string containing the name of the sample where the reads
    are coming from.
  }
}

\details{
  TODO
}

\section{About read names}{
  The read names are typically imported from the BAM file by calling
  \code{\link[GenomicAlignments]{readGAlignments}} (or
  \code{\link[GenomicAlignments]{readGAlignmentPairs}}) with
  \code{use.names=TRUE}. This extracts the "query names" from the
  file (stored in the QNAME field), and makes them the names of the
  returned object.

  The \code{reads} object must have unique names on it. The presence of
  duplicated names generally indicates one (or both) of the following
  situations:

  \itemize{
    \item (a) \code{reads} contains paired-end reads that have not been
              paired;
    \item (b) some of the reads are \emph{secondary alignments}.
  }

  If (a): you can find out whether reads in a BAM file are single- or
  paired-end with the \code{\link[Rsamtools]{quickBamFlagSummary}} utility
  from the \pkg{Rsamtools} package. If they're paired-end, load them
  with \code{\link[GenomicAlignments]{readGAlignmentPairs}}
  instead of \code{\link[GenomicAlignments]{readGAlignments}}, and that
  will pair them.

  If (b): you can filter out secondary alignments by passing
  \code{'isSecondaryAlignment=FALSE'} to \code{\link[Rsamtools]{scanBamFlag}}
  when preparing the \link[Rsamtools]{ScanBamParam} object used to load
  the reads. For example:
  \preformatted{
    flag0 <- scanBamFlag(isSecondaryAlignment=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0)
    reads <- readGAlignments("path/to/BAM/file", use.names=TRUE,
                                  param=param0)
  }
  This will filter out records that have flag 0x100 (secondary alignment)
  set to 1. See \code{?\link[Rsamtools]{scanBamFlag}} in the \pkg{Rsamtools}
  package for more information.
  See the SAM Specs on the SAMtools project page at
  \url{http://samtools.sourceforge.net/} for a description of the
  SAM/BAM flags.
}

\value{
  For \code{assignReads}: the supplied \link{SplicingGraphs} object with
  the reads assigned to it.

  For \code{removeReads}: the supplied \link{SplicingGraphs} object with
  all reads previously assigned with \code{assignReads} removed from it.
}

\author{
  H. Pagès
}

\seealso{
  This man page is part of the \pkg{SplicingGraphs} package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.

  Other topics related to this man page and documented in other packages:
  \itemize{
    \item The \link[GenomicRanges]{GRangesList} class defined in the
          \pkg{GenomicRanges} package.

    \item The \link[GenomicAlignments]{GAlignments} and
          \link[GenomicAlignments]{GAlignmentPairs} classes, as well as
          the \code{\link[GenomicAlignments]{readGAlignments}} and
          \code{\link[GenomicAlignments]{readGAlignmentPairs}} functions
          defined in the \pkg{GenomicAlignments} package.

    \item The \code{\link[Rsamtools]{quickBamFlagSummary}} and
          \code{\link[Rsamtools]{ScanBamParam}} functions defined in the
          \pkg{Rsamtools} package.
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
## 2. Load toy reads
## ---------------------------------------------------------------------
## Load toy reads (single-end) from a BAM file. We filter out secondary
## alignments, reads not passing quality controls, and PCR or optical
## duplicates (see ?scanBamFlag in the Rsamtools package for more
## information):
flag0 <- scanBamFlag(isSecondaryAlignment=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isDuplicate=FALSE)
param0 <- ScanBamParam(flag=flag0)
gal <- readGAlignments(toy_reads_bam(), use.names=TRUE, param=param0)
gal

## ---------------------------------------------------------------------
## 3. Assign the reads to the exons and introns in 'sg'
## ---------------------------------------------------------------------
## The same read can be assigned to more than 1 exon or intron (e.g. a
## junction read with 1 gap can be assigned to 2 exons and 1 intron).
sg <- assignReads(sg, gal, sample.name="TOYREADS")

## See the assignments to the splicing graph edges.
edge_by_tx <- sgedgesByTranscript(sg, with.hits.mcols=TRUE)
edge_data <- mcols(unlist(edge_by_tx))
colnames(edge_data)
head(edge_data)
edge_data[ , c("sgedge_id", "TOYREADS.hits")]

edge_by_gene <- sgedgesByGene(sg, with.hits.mcols=TRUE)
mcols(unlist(edge_by_gene))

## See the assignments to the reduced splicing graph edges.
redge_by_gene <- rsgedgesByGene(sg, with.hits.mcols=TRUE)
mcols(unlist(redge_by_gene))

## ---------------------------------------------------------------------
## 4. Summarize the reads assigned to 'sg' and eventually remove them
## ---------------------------------------------------------------------
## See '?countReads'.
}
