\name{countReads-methods}

\alias{countReads-methods}

\alias{countReads}
\alias{countReads,SplicingGraphs-method}


\title{
  Summarize the reads assigned to a SplicingGraphs object
}

\description{
  \code{getReads} returns the reads assigned to a SplicingGraphs object,
  summarized either by splicing graph edge, \emph{reduced} splicing graph
  edge, transcript, or gene.

  \code{countReads} counts the reads assigned to a SplicingGraphs object.
  The counting can be done by splicing graph edge, \emph{reduced} splicing
  graph edge, transcript, or gene.
}

\usage{
getReads(x, by=c("sgedge", "rsgedge", "tx", "gene"))
countReads(x, by=c("sgedge", "rsgedge", "tx", "gene"))
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object.
  }
  \item{by}{
    Summarize/count by splicing graph edge (\code{by="sgedge"}), by
    \emph{reduced} splicing graph edge (\code{by="rsgedge"}), by transcript
    (\code{by="tx"}), or by gene (\code{by="gene"}).
  }
}

\details{
  TODO
}

\value{
  A \link[IRanges]{DataFrame} object with one row per:
  \itemize{
    \item unique splicing graph edge, if \code{by="sgedge"};
    \item unique \emph{reduced} splicing graph edge, if \code{by="rsgedge"};
    \item transcript if \code{by="tx"};
    \item gene if \code{by="gene"}.
  }

  And with one column per sample (containing the reads for that sample for
  \code{getReads}, and the counts for that sample for \code{countReads}),
  plus the following two additional leading columns:
  \itemize{
    \item if \code{by="sgedge"}: \code{"sgedge_id"}, containing the
          \emph{global splicing graph edge ids}, and \code{"ex_or_in"},
          containing the type of edge (exon or intron);
    \item if \code{by="rsgedge"}: \code{"rsgedge_id"}, containing the
          \emph{global reduced splicing graph edge ids}, and
          \code{"ex_or_in"}, containing the type of edge (exon, intron,
          or mixed);
    \item if \code{by="tx"}: \code{"tx_id"} and \code{"gene_id"};
    \item if \code{by="gene"}: \code{"gene_id"} and \code{"tx_id"}.
  }
}

\author{
  H. Pages
}

\seealso{
  This man page is part of the SplicingGraphs package.
  Please see \code{?`\link{SplicingGraphs-package}`} for an overview of the
  package and for an index of its man pages.
}

\examples{
## ---------------------------------------------------------------------
## 1. Make SplicingGraphs object 'sg' from toy gene model and assign toy
##    reads to it (see '?assignReads')
## ---------------------------------------------------------------------
example(assignReads)

## ---------------------------------------------------------------------
## 2. Summarize the reads by splicing graph edge
## ---------------------------------------------------------------------
getReads(sg) 
countReads(sg)

## ---------------------------------------------------------------------
## 3. Summarize the reads by reduced splicing graph edge
## ---------------------------------------------------------------------
getReads(sg, by="rsgedge")
countReads(sg, by="rsgedge")

## ---------------------------------------------------------------------
## 4. Summarize the reads by transcript
## ---------------------------------------------------------------------
getReads(sg, by="tx")
countReads(sg, by="tx")

## ---------------------------------------------------------------------
## 4. Summarize the reads by gene
## ---------------------------------------------------------------------
getReads(sg, by="gene")
countReads(sg, by="gene")

## ---------------------------------------------------------------------
## 5. Remove the reads from 'sg'.
## ---------------------------------------------------------------------
sg <- removeReads(sg)
countReads(sg)
}
