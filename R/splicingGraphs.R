### =========================================================================
### splicingGraphs() and related utilities
### -------------------------------------------------------------------------

setOldClass("igraph")

.EX_OR_IN_LEVELS2 <- c("ex", "in", "", "mixed")
.EDGE_WEIGHTS <- c(1, 0.2, 0.1, 0.4)
.EX_OR_IN_LEVELS <- .EX_OR_IN_LEVELS2[-4L]


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Spath() accessor
###
### Gets all the splicing paths for the specified gene.
### Returns them in a named IntegerList with 1 top-level element per
### transcript in the specified gene. Each top-level element 'Spath[[i]]'
### contains the splicing site ids for the i-th transcript.
###

setGeneric("Spath", signature="x",
    function(x, gene_id=NA) standardGeneric("Spath")
)

setMethod("Spath", "GRangesList",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        ans <- mcols(x)[ , "Spath"]  # CompressedIntegerList
        if (is.null(x_names)) {
            if (!is.na(gene_id))
                stop("the 'gene_id' arg is not supported ",
                     "when 'x' is unnamed (in which case all its elements ",
                     "(i.e. transcripts) are considered to belong to the ",
                     "same gene)")
            return(ans) 
        }
        if (is.na(gene_id))
            stop("'gene_id' must be supplied when 'x' has names")
        ans <- ans[x_names == gene_id]
        if (length(ans) == 0L)
            stop("invalid 'gene_id'")
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UAhc() accessor
###

setGeneric("UAhc", signature="x",
    function(x, gene_id=NA) standardGeneric("UAhc")
)

setMethod("UAhc", "GRangesList",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        ans <- mcols(x)[["UAhc"]]  # integer vector
        if (is.null(x_names)) {
            if (!is.na(gene_id))
                stop("the 'gene_id' arg is not supported ",
                     "when 'x' is unnamed (in which case all its elements ",
                     "(i.e. transcripts) are considered to belong to the ",
                     "same gene)")
            return(ans) 
        }
        if (is.na(gene_id))
            stop("'gene_id' must be supplied when 'x' has names")
        if (is.null(ans))
            return(ans)
        ans <- ans[x_names == gene_id]
        if (length(ans) == 0L)
            stop("invalid 'gene_id'")
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sgdf() accessor
###
### Returns the splicing graph in a DataFrame with 1 row per edge.
###

### 'spath' must be an IntegerList containing all the splicing paths for a
### given gene. Should have been obtained thru the Spath() accessor.
### Returns a 4-col (or 5-col if 'UAhc' is supplied) data.frame representing
### the splicing graph.
.make_Sgdf0_from_Spath <- function(spath, UAhc=NULL)
{
    if (!is.null(UAhc)) {
        if (!is.integer(UAhc))
            stop("'UAhc' must be an integer vector or NULL")
        if (length(UAhc) != length(spath))
            stop("when not NULL, 'UAhc' must have the same length as 'spath'")
    }
    sgdf0s <- lapply(seq_along(spath),
                     function(i) {
                         SSids <- spath[[i]]
                         from <- c("R", SSids)
                         to <- c(SSids, "L")
                         nb_SSids <- length(SSids)
                         if (nb_SSids %% 2L != 0L)
                             stop("some splicing paths in 'spath' go thru an ",
                                  "odd number of splicing site ids")
                         nexons <- nb_SSids %/% 2L
                         if (nexons == 0L) {
                             ex_or_in <- .EX_OR_IN_LEVELS[3L]
                         } else {
                             nintrons <- nexons - 1L
                             ex_or_in <- c(.EX_OR_IN_LEVELS[3L],
                                           rep.int(.EX_OR_IN_LEVELS[1:2],
                                                   nintrons),
                                           .EX_OR_IN_LEVELS[1L],
                                           .EX_OR_IN_LEVELS[3L])
                         }
                         ex_or_in <- factor(ex_or_in,
                                            levels=.EX_OR_IN_LEVELS)
                         data.frame(from=from,
                                    to=to,
                                    ex_or_in=ex_or_in,
                                    stringsAsFactors=FALSE)
                     })
    nedges_per_tx <- sapply(sgdf0s, nrow)
    sgdf0 <- do.call(rbind, sgdf0s)
    tx_id <- names(spath)
    if (is.null(tx_id))
        tx_id <- seq_along(spath)
    tx_id <- rep.int(factor(tx_id, levels=tx_id), nedges_per_tx)
    sgdf0$tx_id <- tx_id
    if (!is.null(UAhc))
        sgdf0$UAhc <- rep.int(UAhc, nedges_per_tx)
    sgdf0
}

### Collapse the duplicated edges in 'sgdf0' into a DataFrame.
### We use a DataFrame instead of a data.frame because we want to store
### the tx_id col in a CompressedFactorList (even though this container
### doesn't formally exist and a CompressedIntegerList is actually what's
### being used).
.make_Sgdf_from_Sgdf0 <- function(sgdf0)
{
    from <- sgdf0[ , "from"]
    to <- sgdf0[ , "to"]
    ex_or_in <- sgdf0[ , "ex_or_in"]
    tx_id <- sgdf0[ , "tx_id"]
    edges <- paste(from, to, sep="~")
    sm <- match(edges, edges)
    if (!all(ex_or_in == ex_or_in[sm]))
        stop("invalid splicing graph")
    sgdf <- DataFrame(sgdf0[sm == seq_along(sm), , drop=FALSE])
    sgdf$tx_id <- splitAsList(tx_id, sm)
    UAhc <- sgdf$UAhc
    if (!is.null(UAhc))
        sgdf$UAhc <- sum(splitAsList(sgdf0$UAhc, sm))
    sgdf
}

setGeneric("Sgdf", signature="x",
    function(x, gene_id=NA, UAhc=NULL, keep.dup.edges=FALSE)
        standardGeneric("Sgdf")
)

setMethod("Sgdf", "ANY",
    function(x, gene_id=NA, UAhc=NULL, keep.dup.edges=FALSE)
    {
        spath <- Spath(x, gene_id=gene_id)
        if (is.null(UAhc))
            UAhc <- UAhc(x, gene_id=gene_id)
        Sgdf(spath, UAhc=UAhc, keep.dup.edges=keep.dup.edges)
    }
)

setMethod("Sgdf", "IntegerList",
    function(x, gene_id=NA, UAhc=NULL, keep.dup.edges=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an IntegerList")
        sgdf0 <- .make_Sgdf0_from_Spath(x, UAhc=UAhc)
        Sgdf(sgdf0, keep.dup.edges=keep.dup.edges)
    }
)

setMethod("Sgdf", "data.frame",
    function(x, gene_id=NA, UAhc=NULL, keep.dup.edges=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!is.null(UAhc))
            stop("the 'UAhc' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        if (keep.dup.edges)
            return(x)  # no-op
        .make_Sgdf_from_Sgdf0(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various attempts at converting a Splicing graph data frame (Sgdf) into
### something plottable.
###

### ABANDONNED!
### grapNEL objects (defined in the graph package) don't allow more than 1
### edge between the same 2 nodes. Bummer :-(
#.make_graphNEL_from_Sgdf <- function(sgdf)
#{
#    library(graph)
#    from <- sgdf[ , "from"]
#    to <- sgdf[ , "to"]
#    ex_or_in <- sgdf[ , "ex_or_in"]
#    ex_or_in_levels <- levels(ex_or_in)
#    if (!identical(ex_or_in_levels, .EX_OR_IN_LEVELS2)
#     && !identical(ex_or_in_levels, .EX_OR_IN_LEVELS))
#        stop("\"ex_or_in\" column has invalid levels")
#    weights <- .EDGE_WEIGHTS[ex_or_in]
#    nodes <- unique(c(from, to))
#    nodes <- sort(as.integer(setdiff(nodes, c("R", "L"))))
#    nodes <- c("R", as.character(nodes), "L")
#    graph_nel <- graph::graphNEL(nodes, edgemode="directed")
#    graph::addEdge(from, to, graph_nel, weights)
#}

### ABANDONNED!
### Attempt to do the plotting thru an Ragraph object (as a way to circumvent
### the limitation that grapNEL objects don't allow more than 1 edge between
### the same 2 nodes). Pb is that, even though Ragraph objects (defined in
### Rgraphviz) do allow more than 1 edge between the same 2 nodes, the standard
### way to set data on the edges of an Ragraph object is to use `edgeData<-`,
### which unfortunately requires the edges to be specified by the 2 nodes they
### connect. Therefore it is impossible to set different data on the various
### edges between the same 2 nodes, which is a problem if we want to set
### different plotting attributes on them like color, line width, or line
### style.
#.make_Ragraph_from_Sgdf <- function(sgdf, gene_id=NA)
#{
#    library(Rgraphviz)
#    graph_nel <- .make_graphNEL_from_Sgdf(sgdf)
#    pNodes <- Rgraphviz::buildNodeList(graph_nel)
#    pEdges <- Rgraphviz::buildEdgeList(graph_nel, "distinct")
#    from <- sgdf[ , "from"]
#    to <- sgdf[ , "to"]
#    edges <- paste(from, to, sep="~")
#    if (!all(edges %in% names(pEdges)))
#        stop("internal error: the list of pEdge objects returned by ",
#             "buildEdgeList() doesn't have the expected names")
#    pEdges <- pEdges[edges]
#    Rgraphviz::agopen(name=gene_id, nodes=pNodes, edges=pEdges,
#                      edgeMode="directed")
#}

### Experimental. Not ready yet!
.layout.Sgraph <- function(graph)
{
    ## Compute the 'x' col.
    vertices <- get.data.frame(graph, what="vertices")
    nodes <- vertices$name
    nnodes <- length(nodes)
    if (!identical(nodes[c(1L, nnodes)], c("R", "L")))
        stop("first and last nodes are expected to be \"R\" and \"L\", ",
             "respectively")
    max_SSid <- graph$max_SSid
    if (is.null(max_SSid)) {
        nxslot <- nnodes
        xslot <- seq_len(nxslot) - 1L
    } else {
        nxslot <- as.integer(max_SSid) + 2L
        xslot <- c(0L, as.integer(nodes[-c(1L, nnodes)]), nxslot - 1L)
    }
    x <- xslot * 2.0/(nxslot-1L) - 1.0

    ## Compute the 'y' col.
    edges <- get.data.frame(graph, what="edges")
    nedges <- nrow(edges)
    nin <- tabulate(factor(edges$to, levels=nodes), nbins=nnodes)
    nout <- tabulate(factor(edges$from, levels=nodes), nbins=nnodes)
    ## A sanity check.
    if (nin[1L] != 0L || nout[nnodes] != 0L)
        stop("\"R\" or \"L\" nodes cannot have incoming or outgoing edges, ",
             "respectively")
    nin_cumsum <- cumsum(nin)
    nout_cumsum <- cumsum(nout)
    ## A sanity check.
    if (nin_cumsum[nnodes] != nedges || nout_cumsum[nnodes] != nedges)
        stop("The sum of all the incoming edges for all nodes should be ",
             "equal to the total nb of edges. Same for the outgoing edges.")
    ## Nb of edges passing by.
    nbypass <- cumsum(c(0L, nout[-nnodes]) - nin)
    set.seed(33L)
    yslot <- sapply(nbypass+1L, sample, 1L)
    y <- yslot * 2.0/(nbypass+2L) - 1.0

    ## Return the layout matrix.
    cbind(x=x, y=y)
}

### 'sgdf' must be a data.frame as returned by:
###     Sgdf( , keep.dup.edges=TRUE)
### or a DataFrame as returned by:
###     Sgdf( , keep.dup.edges=FALSE)
### Valid extra cols are: "label", "label.color", "lty", "color", "width"
### and "UAhc". They are used to set graphical parameters on the edges.
.precook_igraph_edges_from_Sgdf <- function(sgdf)
{
    required_colnames <- c("from", "to", "ex_or_in", "tx_id")
    extra_colnames <- c("label", "label.color", "lty", "color",
                        "width", "UAhc")
    extract_colnames <- c(required_colnames,
                          intersect(extra_colnames, colnames(sgdf)))
    ans <- sgdf[ , extract_colnames, drop=FALSE]
    ex_or_in <- ans[ , "ex_or_in"]
    ex_or_in_levels <- levels(ex_or_in)
    if (!identical(ex_or_in_levels, .EX_OR_IN_LEVELS2)
     && !identical(ex_or_in_levels, .EX_OR_IN_LEVELS))
        stop("\"ex_or_in\" column has invalid levels")
    if (!("label.color" %in% extract_colnames))
        ans$label.color <- "blue"
    if (!("lty" %in% extract_colnames))
        ans$lty <- c("solid", "solid", "dashed", "solid")[ex_or_in]
    if (!("color" %in% extract_colnames))
        ans$color <- c("green3", "darkgrey", "grey", "black")[ex_or_in]
    if (!("width" %in% extract_colnames) && "UAhc" %in% extract_colnames) {
        min_UAhc <- min(ans$UAhc)
        if (min_UAhc < 0L) {
            warning("'UAhc' column contains negative values. Cannot use ",
                    "it to set the widths of the edges.")
        } else {
            max_UAhc <- max(ans$UAhc)
            if (max_UAhc <= 0L) {
                warning("'UAhc' column has no positive values. Cannot use ",
                        "it to set the widths of the edges.")
            } else {
                ans$width <- 20.0 * ans$UAhc / max(ans$UAhc)
            }
        }
    }
    ans
}

.make_igraph <- function(d)
{
    ## Prepare the 'vertices' argument to pass to graph.data.frame().
    from <- d[ , "from"]
    to <- d[ , "to"]
    nodes <- unique(c(from, to))
    nodes <- sort(as.integer(setdiff(nodes, c("R", "L"))))
    nodes <- c("R", as.character(nodes), "L")
    color <- c("gray", rep.int("white", length(nodes)-2L), "gray")
    label.color <- "black"
    vertices <- data.frame(name=nodes, color=color, label.color=label.color)

    ## Make the igraph object.
    g <- graph.data.frame(d, vertices=vertices)
    layout.kamada.kawai.deterministic <- function(...)
    {
        set.seed(33L)
        layout.kamada.kawai(...)
    }

    ## Set its layout attribute.
    g$layout <- layout.kamada.kawai.deterministic
    #g$layout <- .layout.Sgraph

    g
}

### 'sgdf0' must be a data.frame as returned by:
###     Sgdf( , keep.dup.edges=TRUE)
.make_igraph_from_Sgdf0 <- function(sgdf0, gene_id=NA)
{
    if (!is.data.frame(sgdf0))
        stop("'sgdf0' must be a data.frame")
    d <- .precook_igraph_edges_from_Sgdf(sgdf0)
    if (!("label" %in% colnames(d)))
        d$label <- d$tx_id
    .make_igraph(d)
}

### 'sgdf' must be a DataFrame as returned by:
###     Sgdf( , keep.dup.edges=FALSE)
### or by:
###     Sgdf2( )
.make_igraph_from_Sgdf <- function(sgdf, gene_id=NA)
{
    if (!is(sgdf, "DataFrame"))
        stop("'sgdf' must be a DataFrame")
    d <- .precook_igraph_edges_from_Sgdf(sgdf)
    if (!("label" %in% colnames(d)))
        d$label <- sapply(d$tx_id, paste, collapse=",")
    d$tx_id <- NULL
    ## Turning 'd' into an ordinary data.frame. (Looks like 'as.data.frame()'
    ## on a DataFrame ignores the 'stringsAsFactors' arg so we use
    ## 'data.frame(as.list())' instead.)
    d <- data.frame(as.list(d), stringsAsFactors=FALSE)
    .make_igraph(d)
}

### FIXME: Node and edge attributes are lost. Then plotting is poor...
.make_graphNEL_from_igraph <- function(igraph)
{
    warning("'as.graphNEL' not fully supported yet (node and edge ",
            "attributes are lost)")
    igraph.to.graphNEL(igraph)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sgraph() accessor
###
### Returns the splicing graph in an igraph object.
###

setGeneric("Sgraph", signature="x",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.graphNEL=FALSE)
        standardGeneric("Sgraph")
)

setMethod("Sgraph", "ANY",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.graphNEL=FALSE)
    {
        sgdf <- Sgdf(x, gene_id=gene_id, keep.dup.edges=keep.dup.edges)
        Sgraph(sgdf, as.graphNEL=as.graphNEL)
    }
)

setMethod("Sgraph", "data.frame",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.graphNEL=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a data.frame")
        igraph <- .make_igraph_from_Sgdf0(x)
        Sgraph(igraph, as.graphNEL=as.graphNEL)
    }
)

setMethod("Sgraph", "DataFrame",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.graphNEL=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a DataFrame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a DataFrame")
        igraph <- .make_igraph_from_Sgdf(x)
        Sgraph(igraph, as.graphNEL=as.graphNEL)
    }
)

setMethod("Sgraph", "igraph",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.graphNEL=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an igraph object")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is an igraph object")
        if (!isTRUEorFALSE(as.graphNEL))
            stop("'as.graphNEL' must be TRUE or FALSE")
        if (!as.graphNEL)
            return(x)  # no-op
        .make_graphNEL_from_igraph(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### uninformativeSSids() accessor
###

setGeneric("uninformativeSSids", signature="x",
    function(x, gene_id=NA) standardGeneric("uninformativeSSids")
)

setMethod("uninformativeSSids", "ANY",
    function(x, gene_id=NA)
    {
        x <- Sgdf(x, gene_id=gene_id)
        uninformativeSSids(x)
    }
)

setMethod("uninformativeSSids", "DataFrame",
    function(x, gene_id=NA)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a DataFrame")
        from <- x[ , "from"]
        to <- x[ , "to"]
        from1_SSids <- setdiff(from, from[duplicated(from)])
        to1_SSids <- setdiff(to, to[duplicated(to)])
        intersect(from1_SSids, to1_SSids)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sgdf2() accessor
###
### Same as Sgdf() except that uninformative nodes (i.e. SSids) are removed.
###

### 'sgdf' must be a DataFrame as returned by:
###     Sgdf( , keep.dup.edges=FALSE)
.remove_uninformative_SSids <- function(sgdf)
{
    ex_or_in <- sgdf[ , "ex_or_in"]
    ex_or_in_levels <- levels(ex_or_in)
    if (!identical(ex_or_in_levels, .EX_OR_IN_LEVELS))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, the \"ex_or_in\" column has invalid ",
             "levels. Could it be that it was obtained by a previous call ",
             "to Sgdf2()?")
    levels(ex_or_in) <- .EX_OR_IN_LEVELS2
    uninformative_SSids <- uninformativeSSids(sgdf)
    if (length(uninformative_SSids) == 0L)
        return(sgdf)
    from <- sgdf[ , "from"]
    to <- sgdf[ , "to"]
    tx_id <- sgdf[ , "tx_id"]
    idx1 <- match(uninformative_SSids, from)
    idx2 <- match(uninformative_SSids, to)
    ## 2 sanity checks.
    if (!identical(unname(tx_id[idx1]), unname(tx_id[idx2])))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, the 2 rows containing a given ",
             "uninformative splicing site id must contain the same tx_id.",
             "Could it be that the \"tx_id\" column was manually altered ",
             "before the data.frame (or DataFrame) was passed to ",
             "Sgdf2()?")
    if (!all(idx1 == idx2 + 1L))
        stop("Malformed input.\n",
             "  In the input data.frame (or DataFrame) representing the ",
             "original splicing graph, each uninformative splicing site ",
             "id must appear in 2 consecutive rows (first in the \"to\" ",
             "column, then in the \"from\" column. Could it be that the ",
             "rows were subsetted before the data.frame (or DataFrame) ",
             "was passed to Sgdf2()?")
    from <- from[-idx1]
    to <- to[-idx2]
    ex_or_in[idx1] <- .EX_OR_IN_LEVELS2[4L]
    ex_or_in <- ex_or_in[-idx2]
    tx_id <- tx_id[-idx1]
    DataFrame(from=from, to=to, ex_or_in=ex_or_in, tx_id=tx_id)
}

Sgdf2 <- function(x, gene_id=NA)
{
    if (!is(x, "DataFrame"))
        x <- Sgdf(x, gene_id=gene_id)
    .remove_uninformative_SSids(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sgraph2() accessor
###
### Same as Sgraph() except that uninformative nodes (i.e. SSids) are removed.
###

Sgraph2 <- function(x, gene_id=NA, as.graphNEL=FALSE)
{
    if (!is(x, "DataFrame"))
        x <- Sgdf2(x, gene_id=gene_id)
    Sgraph(x, as.graphNEL=as.graphNEL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### plotSgraph() and plotSgraph2()
###

### ABANDONNED!
### 'x' should be the graphNEL returned by .make_graphNEL_from_Sgdf()
#.plot_graphNEL <- function(x, gene_id=NA)
#{
#    library(Rgraphviz)
#    parameters <- list(nodes=list(col="grey", lty=1, lwd=2,
#                       fontsize=12, height=50,
#                       width=50, fill="#33CC66", shape="ellipse",
#                       fixedsize=FALSE),
#                       graph=list(main=
#                       paste("Splicing graph for gene id:", gene_id)),
#                       edges=list(col="grey", lwd=1.7, labels=NA,
#                       lty=1))
#    x <- Rgraphviz::layoutGraph(x)
#    graph::nodeRenderInfo(x) <- list(fill=c(R="#FF6666", L="#FF6666"))
#    Rgraphviz::renderGraph(x, recipEdges="distinct", graph.pars=parameters,
#                           name="")
#}

plotSgraph <- function(x, gene_id=NA, keep.dup.edges=FALSE)
{
    if (!is(x, "igraph"))
        x <- Sgraph(x, gene_id=gene_id, keep.dup.edges=keep.dup.edges)
    plot.igraph(x, main=gene_id)
}

plotSgraph2 <- function(x, gene_id=NA)
{
    if (!is(x, "igraph"))
        x <- Sgraph2(x, gene_id=gene_id)
    plot.igraph(x, main=gene_id)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### splicingGraphs()
###

### 'exons_start' and 'exons_end' must be 2 integer vectors of the same length
### N representing the start and end positions of N exons that belong to the
### same gene. Returns 2 integer vectors of length N containing the splicing
### site ids assigned to the start and end positions, respectively.
.make_SSids <- function(exons_start, exons_end, on.minus.strand=FALSE)
{
    if (!is.numeric(exons_start))
        stop("'exons_start' must be an integer vector")
    if (!is.integer(exons_start))
        exons_start <- as.integer(exons_start)
    if (!is.numeric(exons_end))
        stop("'exons_end' must be an integer vector")
    if (!is.integer(exons_end))
        exons_end <- as.integer(exons_end)
    N <- length(exons_start)
    if (length(exons_end) != N)
        stop("'exons_start' and 'exons_end' must have the same length")
    ## splicing sites
    ssites <- data.frame(pos=c(exons_start, exons_end),
                         type=rep.int(1:2, c(N, N)))
    ssites_sm <- IRanges:::selfmatchIntegerPairs(ssites$pos, ssites$type)
    ## unique splicing sites
    ussites <- ssites[ssites_sm == seq_along(ssites_sm), , drop=FALSE]
    oo <- IRanges:::orderIntegerPairs(ussites$pos, ussites$type,
                                      decreasing=on.minus.strand)
    ussites <- ussites[oo, , drop=FALSE]
    SSid <- IRanges:::matchIntegerPairs(ssites$pos, ssites$type,
                                        ussites$pos, ussites$type)
    exons_start_SSid <- head(SSid, n=N)
    exons_end_SSid <- tail(SSid, n=N)
    list(start_SSid=exons_start_SSid, end_SSid=exons_end_SSid)
}

### 'exbytx' must be a GRangesList object containing the exons of a *single*
### gene grouped by transcripts. More precisely, each top-level element
### in 'exbytx' contains the genomic ranges of the exons for a particular
### transcript of the gene.
.setSplicingGraphInfo <- function(exbytx, check.introns=TRUE)
{
    if (!is(exbytx, "GRangesList"))
        stop("'exbytx' must be a GRangesList object")
    if (!isTRUEorFALSE(check.introns))
        stop("'check.introns' must be TRUE or FALSE")
    exons <- exbytx@unlistData
    exons_strand <- strand(exons)
    if (nrun(seqnames(exons)) != 1L || nrun(exons_strand) != 1L)
        stop("all the exons in the gene must be on the same ",
             "reference sequence and strand")
    on.minus.strand <- runValue(exons_strand)[1L] == "-"
    if (check.introns) {
        ## We check that, within each transcript, exons are ordered from 5'
        ## to 3' with gaps of at least 1 nucleotide between them.
        ranges_by_tx <- ranges(exbytx)
        if (on.minus.strand)
            ranges_by_tx <- revElements(ranges_by_tx)
        if (!all(isNormal(ranges_by_tx)))
            stop("some transcripts in the gene don't have their exons ",
                 "ordered from 5' to 3' with gaps of at least 1 nucleotide ",
                 "between them (introns)")
    }

    ## Set splicing site ids.
    SSids <- .make_SSids(start(exons), end(exons),
                         on.minus.strand=on.minus.strand)
    exons_mcols <- mcols(exons)
    if (any(names(SSids) %in% colnames(exons_mcols))) {
        in_1_string <- paste(names(SSids), collapse=", ")
        stop("'unlist(exbytx)' already has metadata columns: ", in_1_string)
    }
    mcols(exons) <- cbind(mcols(exons), DataFrame(SSids))
    exbytx@unlistData <- exons
    exbytx_mcols <- mcols(exbytx)

    ## Set tx_id metadata col.
    if ("tx_id" %in% colnames(exbytx_mcols))
        stop("'exbytx' already has metadata column tx_id")
    tx_id <- names(exbytx)
    if (!is.null(tx_id))
        exbytx_mcols$tx_id <- tx_id

    ## Set Spath metadata col.
    if ("Spath" %in% colnames(exbytx_mcols))
        stop("'exbytx' already has metadata column Spath")
    if (on.minus.strand) {
        Spath <- rbind(SSids$end_SSid, SSids$start_SSid)
    } else {
        Spath <- rbind(SSids$start_SSid, SSids$end_SSid)
    }
    Spath_partitioning <- PartitioningByEnd(end(PartitioningByEnd(exbytx)) * 2L)
    names(Spath_partitioning) <- tx_id
    Spath <- splitAsList(as.vector(Spath), Spath_partitioning)
    exbytx_mcols$Spath <- Spath

    mcols(exbytx) <- exbytx_mcols
    exbytx
}

### 'exbytx' must be a GRangesList object containing the exons of one or more
### genes grouped by transcripts. More precisely, each top-level element
### in 'exbytx' contains the genomic ranges of the exons for a particular
### transcript. Typically 'exbytx' will be obtained from a TranscriptDb object
### 'txdb' with 'exonsBy(txdb, by="tx")'.
### 'grouping' is an optional object that represents the grouping by gene of
### the top-level elements (i.e. transcripts) in 'exbytx'. It can be either:
###   (a) Missing (i.e. NULL). In that case, all the transcripts in 'exbytx'
###       are considered to belong to the same gene and the GRangesList object
###       returned by splicingGraphs() will be unnamed.
###   (b) A list of integer or character vectors, or an IntegerList, or a
###       CharacterList object, of length the number of genes to process,
###       and where 'grouping[[i]]' is a valid subscript in 'exbytx' pointing
###       to all the transcripts of the i-th gene.
###   (c) A factor, character vector, or integer vector, of length 'exbytx'
###       with 1 level per gene.
###   (d) A named GRangesList object containing transcripts grouped by genes
###       i.e. each top-level element in 'grouping' contains the genomic ranges
###       of the transcripts for a particular gene. In that case, the grouping
###       is inferred from the tx_id (or alternatively tx_name) metadata
###       column of 'unlist(grouping)' and all the values in that column must
###       be in 'names(exbytx)'.
###       If 'exbytx' was obtained with 'exonsBy(txdb, by="tx")', then the
###       GRangesList object used for grouping would typically be obtained with
###       'transcriptsBy(txdb, by="gene")'.
###   (e) A data.frame or DataFrame with 2 character vector columns: a
###       gene_id column (factor, character vector, or integer vector),
###       and a tx_id (or alternatively tx_name) column. In that case, 'exbytx'
###       must be named and all the values in the tx_id (or tx_name) column
###       must be in 'names(exbytx)'.

.checkOrMakeUniqueGroupingNames <- function(grouping)
{
    grouping_names <- names(grouping)
    if (is.null(grouping_names)) {
        warning("set names on 'grouping' (with 'names(grouping) <- ",
                "seq_along(grouping)') as artificial gene ids")
        names(grouping) <- seq_along(grouping)
        return(grouping)
    }
    if (anyDuplicated(grouping_names))
        stop("'grouping' has duplicated names")
    if (any(c(NA, "") %in% grouping_names))
        stop("'grouping' names contains NAs or \"\" (empty string)")
    grouping
}

### Returns a list or IntegerList or CharacterList. Always *named* with the
### gene ids.
.normargGrouping <- function(grouping, exbytx)
{
    ## (b)
    if (is.list(grouping) || is(grouping, "IntegerList")
     || is(grouping, "CharacterList")) {
        return(.checkOrMakeUniqueGroupingNames(grouping))
    }
    ## (c)
    if (is.factor(grouping) || is.character(grouping) || is.integer(grouping)) {
        if (length(grouping) != length(exbytx))
            stop("when 'grouping' is a factor, character vector, or integer ",
                 "vector, it must have the same length as 'exbytx'")
        return(split(seq_along(exbytx), grouping))
    }
    exbytx_names <- names(exbytx)
    ## (d)
    if (is(grouping, "GRangesList")) {
        if (is.null(exbytx_names))
            stop("when 'grouping' is a GRangesList, 'exbytx' must be named ",
                 "with transcript ids (tx_id) or transcript names (tx_name)")
        grouping <- .checkOrMakeUniqueGroupingNames(grouping)
        mcols <- mcols(grouping@unlistData)
        mcolnames <- colnames(mcols)
        for (colname in c("tx_id", "tx_name")) {
            idx <- which(mcolnames == colname)
            if (length(idx) == 0L)
                next
            if (length(idx) >= 2L)
                stop("'unlist(grouping)' has more than 1 ",
                     colname, " metadata column")
            m <- match(mcols[[idx]], exbytx_names)
            if (any(is.na(m)))
                next
            return(splitAsList(m, PartitioningByEnd(grouping)))
        }
        stop("'unlist(grouping)' has no tx_id or tx_name column, ",
             "or they contain values that are not in 'names(exbytx)'")
    }
    ## (e)
    if (is.data.frame(grouping) || is(grouping, "DataFrame")) {
        if (is.null(exbytx_names))
            stop("when 'grouping' is a data.frame or a DataFrame, 'exbytx' ",
                 "must be named with transcript ids (tx_id) or transcript ",
                 "names (tx_name)")
        grouping_colnames <- colnames(grouping)
        idx <- which(grouping_colnames == "gene_id")
        if (length(idx) != 1L)
            stop("when 'grouping' is a data.frame or a DataFrame, it must ",
                 "have exactly 1 gene_id column")
        gene_id <- grouping[[idx]]
        if (!is.factor(gene_id) && !is.character(gene_id)
         && !is.integer(gene_id))
            stop("'grouping$gene_id' must be a factor, character vector, ",
                 "or integer vector")
        for (colname in c("tx_id", "tx_name")) {
            idx <- which(grouping_colnames == colname)
            if (length(idx) == 0L)
                next
            if (length(idx) >= 2L)
                stop("'grouping' has more than 1 ", colname, " column")
            m <- match(grouping[[idx]], exbytx_names)
            if (any(is.na(m)))
                next
            return(split(m, gene_id))
        }
        stop("'grouping' has no tx_id or tx_name column, ",
             "or they contain values that are not in 'names(exbytx)'")
    }
    stop("invalid 'grouping'")

}

### TODO: Improve handling of invalid genes i.e. provide more details about
### which genes were considered invalid and why.
splicingGraphs <- function(exbytx, grouping=NULL, check.introns=TRUE)
{
    if (!is(exbytx, "GRangesList"))
        stop("'exbytx' must be a GRangesList object")
    if (is.null(grouping)) {
        ans <- .setSplicingGraphInfo(exbytx, check.introns=check.introns)
        names(ans) <- NULL
        return(ans)
    }
    grouping <- .normargGrouping(grouping, exbytx)
    ans <- lapply(seq_along(grouping),
                  function(i) {
                      ii <- grouping[[i]]
                      gene <- exbytx[ii]
                      gene2 <- try(.setSplicingGraphInfo(gene,
                                       check.introns=check.introns),
                                   silent=TRUE)
                      if (inherits(gene2, "try-error"))
                          return(NULL)
                      gene2
                  })
    invalid_genes_idx <- which(sapply(ans, is.null))
    nb_invalid_genes <- length(invalid_genes_idx)
    if (nb_invalid_genes != 0L) {
        warning(nb_invalid_genes, " invalid ",
                ifelse(nb_invalid_genes == 1L, "gene was", "genes were"),
                " dropped")
        grouping <- grouping[-invalid_genes_idx]
        ans <- ans[-invalid_genes_idx]
    }
    ans <- do.call(c, ans)
    grouping_names <- names(grouping)
    if (!is.null(grouping_names)) {
        ans_names <- rep.int(grouping_names, elementLengths(grouping))
        names(ans) <- ans_names
    }
    ans
}

