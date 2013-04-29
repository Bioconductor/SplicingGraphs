\name{countReads-methods}

\alias{countReads-methods}

\alias{countReads}
\alias{countReads,SplicingGraphs-method}


\title{
  Summarize the reads assigned to a SplicingGraphs object
}

\description{
  \code{countReads} counts the reads assigned to a SplicingGraphs object.
  The counting can be done by splicing graph edge, \emph{reduced} splicing
  graph edge, transcript, or gene.

  \code{reportReads} is similar to \code{countReads} but returns right before
  the final counting step, that is, the returned DataFrame contains the reads
  instead of their counts.
}

\usage{
countReads(x, by=c("sgedge", "rsgedge", "tx", "gene"))
reportReads(x, by=c("sgedge", "rsgedge", "tx", "gene"))
}

\arguments{
  \item{x}{
    A \link{SplicingGraphs} object.
  }
  \item{by}{
    Can be \code{"sgedge"}, \code{"rsgedge"}, \code{"tx"}, or \code{"gene"}.
    Specifies the \emph{level of resolution} that summarization should be
    performed at. See Details section below.
  }
}

\details{
  \subsection{Levels of resolution}{
    \code{countReads} and \code{reportReads} allow summarization of the reads
    at different levels of resolution. The level of resolution is determined
    by the type of feature that one chooses via the \code{by} argument.
    The supported resolutions are (from highest to lowest resolution):
    \enumerate{
      \item \code{by="sgedge"} for summarization at the splicing graph edge
            level (i.e. at the exons/intron level);
      \item \code{by="rsgedge"} for summarization at the \emph{reduced}
            splicing graph edge level;
      \item \code{by="tx"} for summarization at the transcript level;
      \item \code{by="gene"} for summarization at the gene level.
    }
  }

  \subsection{Relationship between levels of resolution}{
    There is a parent-child relationship between the features
    corresponding to a given level of resolution (the parent features)
    and those corresponding to a higher level of resolution (the child
    features).

    For example, in the case of the 2 first levels of resolution listed
    above, the parent-child relationship is the following: the parent
    features are the \emph{reduced} splicing graph edges, the child features
    are the splicing graph edges, and each parent feature is obtained by
    merging one or more child features together.
    Similarly, transcripts can be seen as parent features of \emph{reduced}
    splicing graph edges, and genes as parent features of transcripts.
    Note that, the rsgedge/sgedge and gene/tx relationships are one-to-many,
    but the tx/rsgedge relationship is many-to-many because a given edge can
    belong to more than one transcript.

    Finally the parent-child relationships between 2 arbitrary levels of
    resolution is defined by combining the relationships between consecutive
    levels. All possible parent-child relationships are summarized in the
    following table:
    \preformatted{
                    | to: sgedge   | to: rsgedge  | to: tx
      --------------+--------------+--------------+------------
      from: rsgedge | one-to-many  |              |            
      from: tx      | many-to-many | many-to-many |            
      from: gene    | one-to-many  | one-to-many  | one-to-many
    }
  }

  \subsection{Multiple hits and ambiguous reads}{
    An important distinction needs to be made between a read that hits a
    given feature multiple times and a read that hits more than one feature.

    If the former, the read is counted/reported only once for that feature.
    For example, when summarizing at the transcript level, a read is
    counted/reported only once for a given transcript, even if that read
    hits more than one splicing graph edge (or \emph{reduced} splicing graph
    edge) associated with that transcript.

    If the latter, the read is said to be \emph{ambiguous}. An ambiguous read
    is currently counted/reported for each feature where it has a hit.
    This is a temporary situation: in the near future the user will be offered
    options to handle ambiguous reads in different ways.
  }

  \subsection{Ambiguous reads and levels of resolution}{
    A read might be ambiguous at one level of resolution but not at the other.
    Also the number of ambiguous reads is typically affected by the level
    of resolution. However, even though higher resolution generally means
    more ambiguous reads, this is only true when the switch from one level
    of resolution to the other implies a parent-child relationship between
    features that is one-to-many.
    So, based on the above table, this is always true, except when
    switching from using \code{by="tx"} to using \code{by="sgedge"} or
    \code{by="rsgedge"}. In those cases, the switch can produce more
    ambiguities but it can also produce less.
  }
}

\value{
  A \link[IRanges]{DataFrame} object with one row per:
  \itemize{
    \item unique splicing graph edge, if \code{by="sgedge"};
    \item unique \emph{reduced} splicing graph edge, if \code{by="rsgedge"};
    \item transcript if \code{by="tx"};
    \item gene if \code{by="gene"}.
  }

  And with one column per sample (containing the counts for that sample for
  \code{countReads}, and the reads for that sample for \code{reportReads}),
  plus the two following left columns:
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

  For \code{countReads}, each column of counts is of type integer and
  is named after the corresponding sample.
  For \code{reportReads}, each column of reads is a CharacterList object
  and its name is the name of the corresponding sample with the
  \code{".hits"} suffix added to it.
  In both cases, the name of the sample is the name that was passed to
  \code{assignReads} when the reads of a given sample were initially
  assigned. See \code{?\link{assignReads}} for more information.
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
countReads(sg)
reportReads(sg) 

## ---------------------------------------------------------------------
## 3. Summarize the reads by reduced splicing graph edge
## ---------------------------------------------------------------------
countReads(sg, by="rsgedge")
reportReads(sg, by="rsgedge")

## ---------------------------------------------------------------------
## 4. Summarize the reads by transcript
## ---------------------------------------------------------------------
countReads(sg, by="tx")
reportReads(sg, by="tx")

## ---------------------------------------------------------------------
## 5. Summarize the reads by gene
## ---------------------------------------------------------------------
countReads(sg, by="gene")
reportReads(sg, by="gene")

## ---------------------------------------------------------------------
## 6. A close look at ambiguous reads
## ---------------------------------------------------------------------
resolutions <- c("sgedge", "rsgedge", "tx", "gene")

reported_reads <- lapply(resolutions,
    function(by) {
        reported_reads <- reportReads(sg, by=by)
        unlist(reported_reads$TOYREADS.hits)
    })

## The set of reported reads is the same at all levels of resolution:
unique_reported_reads <- lapply(reported_reads, unique)
stopifnot(identical(unique_reported_reads,
                    rep(unique_reported_reads[1], 4)))

## Extract ambigous reads for each level of resolution:
ambiguous_reads <- lapply(reported_reads,
                          function(x) unique(x[duplicated(x)]))
names(ambiguous_reads) <- resolutions
ambiguous_reads

## Reads that are ambiguous at the "rsgedge" level must also be
## ambiguous at the "sgedge" level:
stopifnot(all(ambiguous_reads$rsgedge %in% ambiguous_reads$sgedge))

## However, there is no reason why reads that are ambiguous at the
## "tx" level should also be ambiguous at the "sgedge" or "rsgedge"
## level!

## ---------------------------------------------------------------------
## 7. Remove the reads from 'sg'.
## ---------------------------------------------------------------------
sg <- removeReads(sg)
countReads(sg)
}
