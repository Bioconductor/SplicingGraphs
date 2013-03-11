library(SplicingGraphs)

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

### Slideshow of the graphs.
slideshow(sg)

#sgraph(sg, gene_id="100507433")
#sgraph(sg, gene_id="10362")  # gene official symbol: HMG20B
#sgraph(sg, gene_id="11202")  # gene official symbol: KLK8 (REVIEWED)
#sgraph(sg, gene_id="112724")  # gene official symbol: RDH13 (REVIEWED)
#sgraph(sg, gene_id="1153")  # gene official symbol: CIRBP
#sgraph(sg, gene_id="126259")  # gene official symbol: TMIGD2
#sgraph(sg, gene_id="147650")  # gene official symbol: LINC00085
#sgraph(sg, gene_id="147965")  # gene official symbol: FAM98C


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Plots for CIB3 gene (Entrez ID 117286)
###
### Gene landing page at NCBI: http://www.ncbi.nlm.nih.gov/gene/?term=117286
###

pdf("CIB3-sg.pdf", width=2, height=6)
plot(sgraph(sg, gene_id="117286"))
dev.off()

pdf("CIB3-sg2.pdf", width=2.5, height=6)
plot(sgraph2(sg, gene_id="117286"))
dev.off()

library(Gviz)
# Plotting to the PDF device produces an incomplete plot (ax_track is missing)!
# Looks like a bug in Gviz.
ax_track <- GenomeAxisTrack()
grl <- sg[["117286"]]
### We create 1 track per transcript.
tx_tracks <- lapply(seq_along(grl),
                    function(i) {
                      tx <- grl[i]
                      AnnotationTrack(tx, name=mcols(tx)$tx_id,
                                      fill="orange", shape="box")
                    })
### Put the transcript tracks in an order that is visually more consistent
### with the splicing graph plot.
tx_tracks <- tx_tracks[c(3L, 1L, 2L, 4L)]

pdf("CIB3-gene.pdf", width=6, height=3)
plotTracks(c(list(ax_track), tx_tracks))
dev.off()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Plots for ZNF813 gene (Entrez ID 126017)
###
### Gene landing page at NCBI: http://www.ncbi.nlm.nih.gov/gene/?term=126017

pdf("ZNF813-sg.pdf", width=4, height=6)
plot(sgraph(sg, gene_id="126017"))
dev.off()

pdf("ZNF813-sg2.pdf", width=5, height=6)
plot(sgraph2(sg, gene_id="126017"))
dev.off()

library(Gviz)
# Plotting to the PDF device produces an incomplete plot (ax_track is missing)!
# Looks like a bug in Gviz.
ax_track <- GenomeAxisTrack()
grl <- sg[["126017"]]
### We create 1 track per transcript.
tx_tracks <- lapply(seq_along(grl),
                    function(i) {
                      tx <- grl[i]
                      AnnotationTrack(tx, name=mcols(tx)$tx_id,
                                      fill="orange", shape="box")
                    })
### Put the transcript tracks in an order that is visually more consistent
### with the splicing graph plot.
tx_tracks <- tx_tracks[2:1]

pdf("ZNF813-gene.pdf", width=6, height=3)
plotTracks(c(list(ax_track), tx_tracks))
dev.off()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Entrez ID 10332 as a complex graph! (10 transcripts, 22 splicing sites)
### graphviz layout engine could do a better job though (in particular edge
### 9-17 doesn't need to cross any other edge).

