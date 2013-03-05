\name{sgdf-methods}

\alias{sgdf-methods}

\alias{spath}
\alias{spath,SplicingGraphs-method}

\alias{UATXHcount}
\alias{UATXHcount,SplicingGraphs-method}

\alias{sgdf}
\alias{sgdf,ANY-method}
\alias{sgdf,IntegerList-method}
\alias{sgdf,data.frame-method}

\alias{uninformativeSSids}
\alias{uninformativeSSids,ANY-method}
\alias{uninformativeSSids,DataFrame-method}

\alias{sgdf2}


\title{
  Extract a splicing graph as a data frame
}

\description{
  Extract the splicing graph for a given gene from a \link{SplicingGraphs}
  object and return it as a \link[IRanges]{DataFrame}.
}

\usage{
sgdf(x, gene_id=NA, UATXHcount=NULL, in_by_tx=NULL, keep.dup.edges=FALSE)
sgdf2(x, gene_id=NA)

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

sgdf(sg, gene_id="geneA")
}
