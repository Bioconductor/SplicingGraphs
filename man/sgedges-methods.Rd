\name{sgedges-methods}

\alias{sgedges-methods}

\alias{spath}
\alias{spath,SplicingGraphs-method}

\alias{UATXHcount}
\alias{UATXHcount,SplicingGraphs-method}

\alias{sgedges}
\alias{sgedges,ANY-method}
\alias{sgedges,IntegerList-method}
\alias{sgedges,data.frame-method}

\alias{sgnodes}
\alias{sgnodes,ANY-method}
\alias{sgnodes,IntegerList-method}
\alias{sgnodes,data.frame-method}
\alias{sgnodes,DataFrame-method}

\alias{uninformativeSSids}
\alias{uninformativeSSids,ANY-method}
\alias{uninformativeSSids,DataFrame-method}

\alias{sgedges2}


\title{
  Extract the edges of a splicing graph
}

\description{
  Extract the edges of the splicing graph of a given gene from a
  \link{SplicingGraphs} object and return it as a \link[IRanges]{DataFrame}.
}

\usage{
sgedges(x, gene_id=NA, UATXHcount=NULL, in_by_tx=NULL, keep.dup.edges=FALSE)
sgnodes(x, gene_id=NA)

sgedges2(x, gene_id=NA)

## Related utilities:

spath(x, gene_id=NA)
UATXHcount(x, gene_id=NA)
uninformativeSSids(x, gene_id=NA)
}

\arguments{
  \item{x}{
    TODO
  }
  \item{gene_id}{
    TODO
  }
  \item{UATXHcount}{
    TODO
  }
  \item{in_by_tx}{
    TODO
  }
  \item{keep.dup.edges}{
    TODO
  }
}

\details{
  TODO
}

\value{
  TODO
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \link{SplicingGraphs} class.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.
  }
}

\examples{
example(SplicingGraphs)  # create SplicingGraphs object 'sg'
sg

## 'sg' has 1 element per transcript, and each transcript is
## assigned a name that is the id of the gene it belongs to. All the
## transcripts belonging to the same gene are guaranteed to be
## consecutive elements in 'sg'.
names(sg)

sgedges(sg, gene_id="geneA")
sgnodes(sg, gene_id="geneA")
}
