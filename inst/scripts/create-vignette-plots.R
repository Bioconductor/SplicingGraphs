library(SplicingGraphs)

library(TxDb.Hsapiens.UCSC.hg19.refGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)
tx_by_gn <- transcriptsBy(txdb, by="gene")

### Keep only "clean" genes (i.e. genes with all transcripts on the same
### chromosome and strand).
tx_by_gn <- tx_by_gn[lengths(runValue(seqnames(tx_by_gn))) == 1L]
tx_by_gn <- tx_by_gn[lengths(runValue(strand(tx_by_gn))) == 1L]

### Keep only genes that are on chr19:
is_on_chr19 <- as.character(unlist(runValue(seqnames(tx_by_gn)))) == "chr19"
tx_by_gn <- tx_by_gn[is_on_chr19]

### Compute the splicing graphs (takes about 10 sec.):
sg <- SplicingGraphs(ex_by_tx, tx_by_gn)

### Slideshow of the graphs.
slideshow(sg)

#sgraph(sg["100996665"])  # official gene symbol: LINC01533
#sgraph(sg["55663"])      # official gene symbol: ZNF446 (used in Fig. 2)
#sgraph(sg["100507433"])  # official gene symbol: ZNF571-AS1
#sgraph(sg["1153"])       # official gene symbol: CIRBP
#sgraph(sg["126259"])     # official gene symbol: TMIGD2
#sgraph(sg["147650"])     # official gene symbol: LINC00085
#sgraph(sg["147965"])     # official gene symbol: FAM98C


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Plots for LINC01533 gene (Entrez ID 100996665, 3 transcripts on + strand)
### Used in vignette Fig. 1
###
### Gene landing page at NCBI: http://www.ncbi.nlm.nih.gov/gene/?term=100996665
###

pdf("LINC01533-sg.pdf", width=2, height=6)
plot(sgraph(sg["100996665"]))
dev.off()

pdf("LINC01533-sg2.pdf", width=2.5, height=6)
plot(rsgraph(sg["100996665"]))
dev.off()

grl <- sg[["100996665"]]
### Put the transcripts in an order that is visually more consistent with
### the splicing graph plot.
grl <- grl[c(3L, 2L, 1L)]
pdf("LINC01533-transcripts.pdf", width=6, height=3)
plotTranscripts(grl)
dev.off()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Plots for ZNF446 gene (Entrez ID 55663, 2 transcripts on + strand)
### Used in vignette Fig. 2
### Gene landing page at NCBI: http://www.ncbi.nlm.nih.gov/gene/?term=55663
###

pdf("ZNF446-sg.pdf", width=2, height=6)
plot(sgraph(sg["55663"]))
dev.off()

pdf("ZNF446-sg2.pdf", width=2.5, height=6)
plot(rsgraph(sg["55663"]))
dev.off()

grl <- sg[["55663"]]
### Put the transcripts in an order that is visually more consistent with
### the splicing graph plot.
grl <- grl[c(2L, 1L)]
pdf("ZNF446-transcripts.pdf", width=6, height=2.5)
plotTranscripts(grl)
dev.off()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Plots for ZNF571-AS1 gene (Entrez ID 100507433, 4 transcripts on + strand)
###
### Gene landing page at NCBI: http://www.ncbi.nlm.nih.gov/gene/?term=100507433
###

pdf("ZNF571-AS1-sg.pdf", width=2, height=6)
plot(sgraph(sg["100507433"]))
dev.off()

pdf("ZNF571-AS1-sg2.pdf", width=2.5, height=6)
plot(rsgraph(sg["100507433"]))
dev.off()

grl <- sg[["100507433"]]
### Put the transcripts in an order that is visually more consistent with
### the splicing graph plot.
grl <- grl[c(4L, 3L, 1L, 2L)]
pdf("ZNF571-AS1-transcripts.pdf", width=6, height=3)
plotTranscripts(grl)
dev.off()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Entrez ID 10332 as a complex graph! (12 transcripts on + strand, 23
### splicing sites)

