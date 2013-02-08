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

    \item \code{\link{splicingGraphs}} and related tools for creating and
          plotting splicing graphs from transcript annotations.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. Load toy reads and a toy gene model, and find the hits that are
##    "compatible" with the splicing of the transcripts
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

## Load a toy gene model as a TranscriptDb object. It is important to
## use a gene model based on the exact same reference genome that was
## used to align the reads.
toy_genes_txdb <- makeTranscriptDbFromGFF(toy_genes_gff())

## Extract the exons grouped by transcript:
ex_by_tx <- exonsBy(toy_genes_txdb, by="tx", use.names=TRUE)

## Find the overlaps between the reads and the transcripts:
grl <- grglist(gal, order.as.in.query=TRUE)
ov0 <- findOverlaps(grl, ex_by_tx, ignore.strand=TRUE)

## Nb of hits per read, per transcript:
head(countQueryHits(ov0))  # per read
head(countSubjectHits(ov0))  # per transcript

## Keep only overlaps that are "compatible" with the splicing of the
## transcripts:
ovenc0 <- encodeOverlaps(grl, ex_by_tx, hits=ov0,
                         flip.query.if.wrong.strand=TRUE)
ov0_is_comp <- isCompatibleWithSplicing(ovenc0)
ov1 <- ov0[ov0_is_comp]

## ---------------------------------------------------------------------
## 2. Count unambiguous hits per transcript (UnAmbiguous hit counts)
## ---------------------------------------------------------------------

ov2 <- ov1[queryHits(ov1) \%in\% which(countQueryHits(ov1) == 1L)]
UAhc <- countSubjectHits(ov2)
mcols(ex_by_tx)$UAhc <- UAhc

## Compute the splicing graphs (1 per gene), and plot them:
tx_by_gn <- transcriptsBy(toy_genes_txdb, by="gene")
ex_by_tx2 <- splicingGraphs(ex_by_tx, tx_by_gn)
names(ex_by_tx2)
plotSgraph(ex_by_tx2, "geneA")
}
