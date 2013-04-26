\name{countReads-methods}

\alias{countReads-methods}

\alias{countReads}
\alias{countReads,SplicingGraphs-method}


\title{
  Summarize the reads assigned to the edges of a SplicingGraphs object
}

\description{
  \code{countReads} returns a summarized count of the reads assigned to
  the edges of a SplicingGraphs object.
}

\usage{
countReads(x, by=c("sgedge", "rsgedge", "tx", "gene"))
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object.
  }
  \item{by}{
    Summarize by splicing graph edge (\code{by="sgedge"}), by \emph{reduced}
    splicing graph edge (\code{by="rsgedge"}), by transcript (\code{by="tx"}),
    or by gene (\code{by="gene"}).
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

  And with one column per sample (containing the counts for that sample),
  plus the two following additional leading columns:
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
## 2. Summarize the reads assigned to 'sg'
## ---------------------------------------------------------------------
countReads(sg)  # nb of reads per splicing graph edge
countReads(sg, by="rsgedge")  # ... per reduced splicing graph edge
countReads(sg, by="tx")  # ... per transcript
countReads(sg, by="gene")  # ... per gene

## ---------------------------------------------------------------------
## 3. Remove the reads from 'sg'.
## ---------------------------------------------------------------------
removeReads(sg)
countReads(sg)
}
