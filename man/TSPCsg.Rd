\name{TSPCsg}

\alias{TSPCsg}
\alias{TSPC}

\title{
  TSPC splicing graphs
}

\description{
  TODO
}

\examples{
## Load SplicingGraphs object 'TSPCsg':
filepath <- system.file("extdata", "TSPCsg.rda", package="SplicingGraphs")
load(filepath)
TSPCsg

## 'TSPCsg' has 1 element per gene and 'names(sg)' gives the gene ids.
names(TSPCsg)

## 1 splicing graph per gene. (Note that gene MUC16 was dropped
## because transcripts T-4 and T-5 in this gene both have their
## 2nd exon *inside* their 3rd exon. Splicing graph theory doesn't
## apply in that case.)

## Extract the edges of a given graph:
TSPCsgedges <- sgedges(TSPCsg["LGSN"])
TSPCsgedges

## Plot the graph for a given gene:
plot(TSPCsg["LGSN"])  # or 'plot(sgraph(TSPCsgedges))'

## The reads from all samples have been assigned to 'TSPCsg'.
## Use countReads() to summarize by splicing graph edge:
counts <- countReads(TSPCsg)
dim(counts)
counts[ , 1:5]

## You can subset 'TSPCsg' by 1 or more gene ids before calling
## countReads() in order to summarize only for those genes:
DAPL1counts <- countReads(TSPCsg["DAPL1"])
dim(DAPL1counts)
DAPL1counts[ , 1:5]

## Use 'by="rsgedge"' to summarize by *reduced* splicing graph edge:
DAPL1counts2 <- countReads(TSPCsg["DAPL1"], by="rsgedge")
dim(DAPL1counts2)
DAPL1counts2[ , 1:5]

## No reads assigned to genes KIAA0319L or TREM2 because no
## BAM files were provided for those genes:
KIAA0319Lcounts <- countReads(TSPCsg["KIAA0319L"])
KIAA0319Lcountsums <- sapply(KIAA0319Lcounts[ , -(1:2)], sum)
stopifnot(all(KIAA0319Lcountsums == 0))

TREM2counts <- countReads(TSPCsg["TREM2"])
TREM2countsums <- sapply(TREM2counts[ , -(1:2)], sum)
stopifnot(all(TREM2countsums == 0))

## Plot all the splicing graphs:
slideshow(TSPCsg)
}
