\name{TSPCsgdfs}

\alias{TSPCsgdfs}
\alias{TSPC}

\title{
  Loading an plotting the TSPC splicing graph data frames
}

\examples{
## 1 splicing graph data frame per gene, except for gene MUC16.
## Transcripts T-4 and T-5 in gene MUC16 both have their 2nd exon included in
## their 3rd exon ==> splicing graph theory doesn't apply.
filepaths <- list.files(system.file("extdata", "TSPCsgdfs",
                                    package="SplicingGraphs"),
                        full.names=TRUE)
for (filepath in filepaths)
    load(filepath)

dim(BAI1sgdf)

## All the data frames have 1 row per edge in the graph, and the first 4
## cols are always "from", "to", "ex_or_in", and "tx_id". Note that there can
## be more than 1 transcript associated with a given edge.
LGSNsgdf[ , 1:4]

## There is 1 additional column per sample:
LGSNsgdf[ , 5:8]

## 'KIAA0319Lsgdf' and 'TREM2sgdf' have no additional cols because there was
## no BAM files for those genes:
dim(KIAA0319Lsgdf)
dim(TREM2sgdf)

## Plot the splicing graphs:
library(Rgraphviz)
plot(Sgraph(BAI1sgdf, as.graphNEL=TRUE))
plot(Sgraph(CYB561sgdf, as.graphNEL=TRUE))
plot(Sgraph(DAPL1sgdf, as.graphNEL=TRUE))
plot(Sgraph(ITGB8sgdf, as.graphNEL=TRUE))
plot(Sgraph(KIAA0319Lsgdf, as.graphNEL=TRUE))
plot(Sgraph(LGSNsgdf, as.graphNEL=TRUE))
plot(Sgraph(MKRN3sgdf, as.graphNEL=TRUE))
plot(Sgraph(ST14sgdf, as.graphNEL=TRUE))
plot(Sgraph(TREM2sgdf, as.graphNEL=TRUE))
}
