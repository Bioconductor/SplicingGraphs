### =========================================================================
### Unit test for R/countReads-methods.R
### -------------------------------------------------------------------------

.make_toy_sg <- function()
{
     library(GenomicFeatures)
     suppressWarnings(
       toy_genes_txdb <- makeTxDbFromGFF(toy_genes_gff())
     )
     SplicingGraphs(toy_genes_txdb)
}

.load_toy_reads <- function()
{
    flag0 <- scanBamFlag(isSecondaryAlignment=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isDuplicate=FALSE)
    param0 <- ScanBamParam(flag=flag0)
    readGAlignments(toy_reads_bam(), use.names=TRUE, param=param0)
}

test_countReads <- function()
{
    sg <- .make_toy_sg()
    reads <- .load_toy_reads()
    sg <- assignReads(sg, reads, sample.name="TOYREADS")
    ## Shallow testing.
    checkIdentical(c("sgedge_id", "ex_or_in", "TOYREADS"),
                   colnames(countReads(sg[0])))
    checkIdentical(c("rsgedge_id", "ex_or_in", "TOYREADS"),
                   colnames(countReads(sg[0], by="rsgedge")))
    checkIdentical(c("tx_id", "gene_id", "TOYREADS"),
                   colnames(countReads(sg[0], by="tx")))
    checkIdentical(c("gene_id", "tx_id", "TOYREADS"),
                   colnames(countReads(sg[0], by="gene")))
    ## TODO: Deep testing.
}

