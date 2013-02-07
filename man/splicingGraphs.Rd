\name{splicingGraphs}

\alias{Spath}
\alias{Spath,GRangesList-method}

\alias{Sgdf}
\alias{Sgdf,ANY-method}
\alias{Sgdf,IntegerList-method}
\alias{Sgdf,data.frame-method}

\alias{Sgraph}
\alias{Sgraph,ANY-method}
\alias{Sgraph,data.frame-method}
\alias{Sgraph,DataFrame-method}
\alias{Sgraph,igraph-method}

\alias{uninformativeSSids}
\alias{uninformativeSSids,ANY-method}
\alias{uninformativeSSids,DataFrame-method}

\alias{Sgdf2}

\alias{Sgraph2}

\alias{plotSgraph}
\alias{plotSgraph2}

\alias{splicingGraphs}

\title{
  Tools for creating and plotting splicing graphs from transcript annotations
}

\description{
  TODO
}

\usage{
splicingGraphs(exbytx, grouping=NULL, check.introns=TRUE)

## Related utilities
Spath(x, gene_id=NA)
Sgdf(x, gene_id=NA, keep.dup.edges=FALSE)
Sgraph(x, gene_id=NA, keep.dup.edges=FALSE, as.graphNEL=FALSE)
uninformativeSSids(x, gene_id=NA)
Sgdf2(x, gene_id=NA)
Sgraph2(x, gene_id=NA, as.graphNEL=FALSE)
plotSgraph(x, gene_id=NA, keep.dup.edges=FALSE)
plotSgraph2(x, gene_id=NA)
}

\arguments{
  \item{exbytx}{
    A \link[GenomicRanges]{GRangesList} object containing the exons of one
    or more genes grouped by transcripts. See Details section below.
  }
  \item{grouping}{
    An optional object that represents the grouping by gene of the top-level
    elements (i.e. the transcripts) in \code{exbytx}. See Details section
    below.
  }
  \item{check.introns}{
    TODO
  }
  \item{x}{
    TODO
  }
  \item{gene_id}{
    TODO
  }
  \item{keep.dup.edges}{
    TODO
  }
  \item{as.graphNEL}{
    TODO
  }
}

\details{
  \code{exbytx} must be a \link[GenomicRanges]{GRangesList} object
  containing the exons of one or more genes grouped by transcripts. More
  precisely, each top-level element in \code{exbytx} contains the genomic
  ranges of the exons for a particular transcript. Typically \code{exbytx}
  will be obtained from a \link[GenomicFeatures]{TranscriptDb} object
  \code{txdb} with \code{\link[GenomicFeatures]{exonsBy}(txdb, by="tx")}.

  \code{grouping} is an optional object that represents the grouping by gene
  of the top-level elements (i.e. the transcripts) in \code{exbytx}. It can
  be either:
  \itemize{
    \item Missing (i.e. \code{NULL}). In that case, all the transcripts in
          \code{exbytx} are considered to belong to the same gene.
    \item A list of integer or character vectors, or an
          \link[IRanges]{IntegerList}, or a \link[IRanges]{CharacterList}
          object, of length the number of genes to process, and where
          \code{grouping[[i]]} is a valid subscript in \code{exbytx}
          pointing to all the transcripts of the i-th gene.
    \item A factor, character vector, or integer vector, of the same length
          as \code{exbytx} and 1 level per gene.
    \item A named \link[GenomicRanges]{GRangesList} object containing
          transcripts grouped by genes i.e. each top-level element in
          \code{grouping} contains the genomic ranges of the transcripts
          for a particular gene. In that case, the grouping is inferred from
          the tx_id (or alternatively tx_name) metadata column of
          \code{unlist(grouping)} and all the values in that column must
          be in \code{names(exbytx)}.
          If \code{exbytx} was obtained with
          \code{\link[GenomicFeatures]{exonsBy}(txdb, by="tx")},
          then the \link[GenomicRanges]{GRangesList} object used for grouping
          would typically be obtained with
          \code{\link[GenomicFeatures]{transcriptsBy}(txdb, by="gene")}.
    \item A data.frame or \link[IRanges]{DataFrame} with 2 character vector
          columns: a gene_id column (factor, character vector, or integer
          vector), and a tx_id (or alternatively tx_name) column. In that
          case, \code{exbytx} must be named and all the values in the tx_id
          (or tx_name) column must be in \code{names(exbytx)}.
  }
}

\value{
TODO
}

\author{
  H. Pages
}

\references{
  Heber, S., Alekseyev, M., Sze, S., Tang, H., and Pevzner, P. A. 
  \emph{Splicing graphs and EST assembly problem}
  Bioinformatics
  Date: Jul 2002
  Vol: 18
  Pages: S181-S188

  Sammeth, M. (2009) 
  \emph{Complete alternative splicing events are bubbles in splicing graphs}
  J. Comput. Biol.
  Date: Aug 2009
  Vol: 16
  Pages: 1117-1140
}

\seealso{
  \itemize{
    \item The \link[GenomicRanges]{GRangesList} class in the GenomicRanges
          package.
    \item \code{\link[GenomicFeatures]{exonsBy}},
          \code{\link[GenomicFeatures]{transcriptsBy}}, and the
          \link[GenomicFeatures]{TranscriptDb} class in the GenomicFeatures
          package.
    \item The \link[IRanges]{IntegerList}, \link[IRanges]{CharacterList},
          and \link[IRanges]{DataFrame} classes in the IRanges package.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. Load a gene model (made of 2 toy genes) as a TranscriptDb object,
##    and plot it
## ---------------------------------------------------------------------

toy_genes_gff <- system.file("extdata", "toy_genes.gff3", package="SpliceGraph")

## Note that you can display the content of the file with
cat(readLines(toy_genes_gff), sep="\n")

suppressWarnings(
  toy_genes_txdb <- makeTranscriptDbFromGFF(toy_genes_gff)
)

## Plot the gene model
if (interactive()) {
  library(Gviz)
  grtrack <- GeneRegionTrack(toy_genes_txdb, name="toy genes")
  gatrack <- GenomeAxisTrack()
  plotTracks(list(gatrack, grtrack))
}

## ---------------------------------------------------------------------
## 2. Compute all the splicing graphs (1 per gene) into a single
##    GRangesList object
## ---------------------------------------------------------------------

## Extract the exons grouped by transcript
ex_by_tx <- exonsBy(toy_genes_txdb, by="tx")

## Extract the transcript grouped by gene
tx_by_gn <- transcriptsBy(toy_genes_txdb, by="gene")

ex_by_tx2 <- splicingGraphs(ex_by_tx, tx_by_gn)
## Note that 'ex_by_tx2' can also be created directly from the
## TranscriptDb object:
#ex_by_tx2 <- splicingGraphs(toy_genes_txdb)  # not ready yet

## 'ex_by_tx2' has 1 element per transcript and its names are the ids
## of the corresponding gene:
names(ex_by_tx2)

## ---------------------------------------------------------------------
## 3. Create the splicing graph data frame for geneA, and plot it
## ---------------------------------------------------------------------

sgdfA <- Sgdf(ex_by_tx2, gene_id="geneA")
sgdfA

if (interactive()) {
  plotSgraph(sgdfA)
}

## ---------------------------------------------------------------------
## 4. Extract information from the splicing graph
## ---------------------------------------------------------------------

if (FALSE) {
  ## FIXME: All the stuff below is broken!
  ## get the new disjoint exon ids
  mcols(unlist(sgA))$disJ_exon_id

  ## get the gene ids
  mcols(sgA)$gene_id

  ## get the original exon ids
  mcols(unlist(sgA))$exon_ids

  ## access the splice events via the metadata slot
  spliceEvents <- metadata(sgA)$spliceEvents
  spliceEvents

  ## count the splice events
  table(spliceEvents$eventType)

  ## assignment of edges to bubbles
  mcols(sgA)
}

## ---------------------------------------------------------------------
## 5. Repeating 3. for geneB
## ---------------------------------------------------------------------

sgdfB <- Sgdf(ex_by_tx2, gene_id="geneB")
sgdfB

if (interactive()) {
  plotSgraph(sgdfB)
}
}
