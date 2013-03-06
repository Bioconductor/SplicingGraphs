\name{TSPCsgraphs}

\alias{TSPCsgraphs}
\alias{TSPC}

\title{
  Loading an plotting the TSPC splicing graph data frames
}

\description{
  TODO
}

\examples{
## 1 splicing graph data frame per gene, except for gene MUC16.
## Transcripts T-4 and T-5 in gene MUC16 both have their 2nd exon included
## in their 3rd exon ==> splicing graph theory doesn't apply.
filepaths <- list.files(system.file("extdata", "TSPCsgraphs",
                                    package="SplicingGraphs"),
                        full.names=TRUE)
for (filepath in filepaths)
    load(filepath)

dim(BAI1sgedges)

## All the data frames have 1 row per edge in the graph, and the first 4
## cols are always "from", "to", "ex_or_in", and "tx_id". Note that there
## can be more than 1 transcript associated with a given edge.
LGSNsgedges[ , 1:4]

## There is 1 additional column per sample:
LGSNsgedges[ , 5:8]

## 'KIAA0319Lsgedges' and 'TREM2sgedges' have no additional cols because there
## was no BAM files for those genes:
dim(KIAA0319Lsgedges)
dim(TREM2sgedges)

## Plot the splicing graphs:
library(Rgraphviz)
plot(sgraph(BAI1sgedges))
plot(sgraph(CYB561sgedges))
plot(sgraph(DAPL1sgedges))
plot(sgraph(ITGB8sgedges))
plot(sgraph(KIAA0319Lsgedges))
plot(sgraph(LGSNsgedges))
plot(sgraph(MKRN3sgedges))
plot(sgraph(ST14sgedges))
plot(sgraph(TREM2sgedges))
}
