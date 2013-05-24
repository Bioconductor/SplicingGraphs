\name{SplicingGraphs-package}
\docType{package}

\alias{SplicingGraphs-package}

\title{
  Create, manipulate, visualize splicing graphs, and assign RNA-seq reads
  to them
}

\description{
  The \pkg{SplicingGraphs} package allows the user to create, manipulate,
  and visualize splicing graphs and their bubbles based on a gene model for
  a given organism. Additionally it allows the user to assign RNA-seq
  reads to the edges of a set of splicing graphs, and to summarize them in
  different ways.
}

\details{
  See the \emph{Splicing graphs and RNA-seq data} vignette in the package
  for a gentle introduction to its use.
  To access the vignette, do \code{browseVignettes("SplicingGraphs")}, then
  click on the link to the PDF version.
}

\note{
  The \pkg{SplicingGraphs} package is a reincarnation of an internal project,
  the \pkg{SpliceGraph} package, originally written by D. Bindreither,
  M. Carlson, and M. Morgan. \pkg{SpliceGraph} was never released as part of
  Bioconductor.

  With respect to the old \pkg{SpliceGraph}, the scope of the new
  \pkg{SplicingGraphs} package has been reduced to focus only on the following
  functionalities: creating/manipulating/plotting splicing graphs, computing
  the bubbles and AS codes, and assigning/counting reads.

  In addition to this, the old \pkg{SpliceGraph} package also had facilities
  for performing some downstream statistical analysis. They were covered in
  its vignette under the following topics/sections:
  \itemize{
    \item Experimental design
    \item Significant altered alternative splice events
    \item Modification of GLM employed in the \pkg{DEXSeq} package
    \item Differential edge expression analysis
    \item Comparison to the classic \pkg{DEXSeq} analysis
  }
  The \pkg{SplicingGraphs} vignette doesn't cover any of this and the new
  package provides no facilities for doing this type of downstream statistical
  analysis. 
}

\author{
  Author and maintainer: H. Pages <hpages@fhcrc.org>

  The \pkg{SplicingGraphs} package is a complete revamp (design and
  implementation) of the old \pkg{SpliceGraph} package (see Note above).
}

\references{
  Heber, S., Alekseyev, M., Sze, S., Tang, H., and Pevzner, P. A. 
  \emph{Splicing graphs and EST assembly problem}
  Bioinformatics
  Date: Jul 2002
  Vol: 18
  Pages: S181-S188

  Sammeth, M. (2009) 
  \emph{Complete Alternative Splicing Events Are Bubbles in Splicing Graphs}
  J. Comput. Biol.
  Date: Aug 2009
  Vol: 16
  Pages: 1117-1140
}

\seealso{
  The man pages in the \pkg{SplicingGraphs} package are:
  \enumerate{
    \item The \link{SplicingGraphs} class.

    \item \code{\link{plotTranscripts}} for plotting a set of transcripts
          along genomic coordinates.

    \item \code{\link{sgedgesByGene}} for extracting the edges and their
          ranges from a SplicingGraphs object.

    \item \code{\link{txpath}} for extracting the transcript paths of a
          splicing graph.

    \item \code{\link{sgedges}} for extracting the edges (and nodes) of a
          splicing graph.

    \item \code{\link{sgraph}} for extracting a splicing graph as a
          plottable graph-like object.

    \item \code{\link{rsgedgesByGene}} for extracting the reduced edges
          and their ranges from a SplicingGraphs object.

    \item \code{\link{bubbles}} for computing the bubbles of a splicing graph.

    \item \code{\link{assignReads}} for assigning reads to the edges of a
          SplicingGraphs object.

    \item \code{\link{countReads}} for summarizing the reads assigned to
          a SplicingGraphs object.

    \item \code{\link{toy_genes_gff}} for details about the toy data included
          in this package.
  }
}

\examples{
if (interactive()) {
  ## Access the "Splicing graphs and RNA-seq data" vignette with:
  browseVignettes("SplicingGraphs")
}
}
\keyword{package}
