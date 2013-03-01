### Create and plot splicing graph for Human gene JUNB (Entrez ID 3726).

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)
tx_by_gn <- transcriptsBy(txdb, by="gene")

### Keep only clean genes (i.e. genes with all transcripts on the same
### chromosome and strand).
tx_by_gn <- tx_by_gn[elementLengths(runValue(seqnames(tx_by_gn))) == 1L]
tx_by_gn <- tx_by_gn[elementLengths(runValue(strand(tx_by_gn))) == 1L]

### Keep only genes that are on chr19:
is_on_chr19 <- as.character(unlist(runValue(seqnames(tx_by_gn)))) == "chr19"
tx_by_gn <- tx_by_gn[is_on_chr19]

### Keep only genes that have at least 2 transcripts:
tx_by_gn <- tx_by_gn[elementLengths(tx_by_gn) >= 2L]

### Compute the splicing graphs (takes about 2 min. 30 sec.):
sg <- SplicingGraphs(ex_by_tx, tx_by_gn)

for (gene_id in unique(names(sg))) {
    ntx <- sum(names(sg) == gene_id)
    cat("Plotting gene ", gene_id, " (", ntx, " transcripts). ", sep="")
    plot(sg, gene_id)
    cat("Press <Enter> for next...")
    readLines(n=1)
}

### Only VALIDATED genes:
#sgraph(sg, gene_id="100507433")
#sgraph(sg, gene_id="10362")  # gene official symbol: HMG20B
#sgraph(sg, gene_id="11202")  # gene official symbol: KLK8 (REVIEWED)
#sgraph(sg, gene_id="112724")  # gene official symbol: RDH13 (REVIEWED)
#sgraph(sg, gene_id="1153")  # gene official symbol: CIRBP

sgCIB3 <- sgraph(sg, gene_id="117286")  # gene official symbol: CIB3 (REVIEWED)
pdf("CIB3-sg.pdf", width=120, height=480)
plot(sgCIB3)
dev.off()

#sgraph(sg, gene_id="126259")  # gene official symbol: TMIGD2
#sgraph(sg, gene_id="147650")  # gene official symbol: LINC00085

#sgFAM98C <- sgraph(sg, gene_id="147965")  # gene official symbol: FAM98C
#pdf("FAM98C-sg.pdf", width=240, height=480)
#plot(sgFAM98C)
#dev.off()

library(Gviz)
# Plotting to the PDF device produces an incomplete plot (ax_track is missing)!
# Looks like a bug in Gviz.
#pdf("CIB3-gene.pdf", width=480, height=180)
ax_track <- GenomeAxisTrack()
track_list <- list(ax_track)
grl <- sg@tx[names(sg@tx) == "117286"]
for (i in seq_along(grl)) {
    tx <- grl[i]
    tx_track <- AnnotationTrack(tx, name=mcols(tx)$tx_id,
                                fill="orange", shape="box")
    track_list <- c(track_list, list(tx_track))
}
plotTracks(track_list)
#dev.off()

