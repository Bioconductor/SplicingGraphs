\name{spliceGraphs}

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

\examples{
## Load toy reads (single-end) from a BAM file. We filter out secondary
## alignments, reads not passing quality controls, and PCR or optical
## duplicates (see ?scanBamFlag for more information):
flag0 <- scanBamFlag(isNotPrimaryRead=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isDuplicate=FALSE)
param0 <- ScanBamParam(flag=flag0)
gal <- readGappedAlignments(toy_reads_bam(), use.names=TRUE, param=param0)
gal

## Load a toy gene model as a TranscriptDb object:
toy_genes_txdb <- makeTranscriptDbFromGFF(toy_genes_gff())

## Extract the exons grouped by transcript:
ex_by_tx <- exonsBy(toy_genes_txdb, by="tx", use.names=TRUE)

## Find the overlaps between the reads and the transcripts:
grl <- grglist(gal, order.as.in.query=TRUE)
ov0 <- findOverlaps(grl, ex_by_tx, ignore.strand=TRUE)

## Nb of hits per read, per transcript:
head(countQueryHits(ov0))  # per read
head(countSubjectHits(ov0))  # per transcript

## Keep only overlaps "compatible" with the splicing of the transcript:
ovenc0 <- encodeOverlaps(grl, ex_by_tx, hits=ov0,
                         flip.query.if.wrong.strand=TRUE)
ov0_is_comp <- isCompatibleWithSplicing(ovenc0)
ov <- ov0[ov0_is_comp]
}
