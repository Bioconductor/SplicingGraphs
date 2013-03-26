\name{countReads}

\alias{countReads}

\title{
  Tools for assigning read counts to splicing graphs
}

\description{
  TODO
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link[IRanges]{findOverlaps}} and
          \code{\link[IRanges]{encodeOverlaps}} in the IRanges
          package for finding overlapping ranges and encode the overlaps.

    \item \code{\link[GenomicFeatures]{exonsBy}},
          \code{\link[GenomicFeatures]{transcriptsBy}}, and the
          \link[GenomicFeatures]{TranscriptDb} class in the GenomicFeatures
          package.

    \item \code{\link{SplicingGraphs}} and related tools for creating and
          plotting splicing graphs from transcript annotations.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. Make SplicingGraphs object 'sg' from toy gene model in
##    toy_genes_gff() (see '?SplicingGraphs')
## ---------------------------------------------------------------------

example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(sg)

## ---------------------------------------------------------------------
## 2. Load toy reads and find the hits that are "compatible" with the
##    splicing of the transcripts
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

## Assign the reads to the exons/introns in 'sg'. The same read can be
## assigned to more than 1 exon or intron (e.g. a junction read with
## 1 gap can be assigned to 2 exons and 1 intron).
sg <- assignReads(sg, grl)

sg[["geneA"]]
unlist(sg)
intronsByTranscript(sg)
sgedges(sg["geneA"])

## ---------------------------------------------------------------------
## 3. Count nb of unambiguous compatible hits per transcript
## ---------------------------------------------------------------------

#ov2 <- ov1[queryHits(ov1) \%in\% which(countQueryHits(ov1) == 1L)]
#UATXHcount <- countSubjectHits(ov2)
#mcols(sg@genes@unlistData)$UATXHcount <- UATXHcount

#plot(sg["geneA"])
}
