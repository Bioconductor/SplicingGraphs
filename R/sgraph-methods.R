### =========================================================================
### sgraph (and related) methods
### -------------------------------------------------------------------------


setOldClass("igraph")

.EDGE_WEIGHTS <- c(1, 0.2, 0.1, 0.4)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_igraph_from_sgedges()
###

### 'sgedges' must be a data.frame as returned by:
###     sgedges( , keep.dup.edges=TRUE)
### or a DataFrame as returned by:
###     sgedges( , keep.dup.edges=FALSE)
### Valid extra cols are: "label", "label.color", "lty", "color", "width"
### and "UATXHcount". They are used to set graphical parameters on the edges.
.precook_igraph_edges_from_sgedges <- function(sgedges)
{
    required_colnames <- c("from", "to", "ex_or_in", "tx_id")
    extra_colnames <- c("label", "label.color", "lty", "color",
                        "width", "UATXHcount")
    extract_colnames <- c(required_colnames,
                          intersect(extra_colnames, colnames(sgedges)))
    ans <- sgedges[ , extract_colnames, drop=FALSE]
    ex_or_in <- ans[ , "ex_or_in"]
    ex_or_in_levels <- levels(ex_or_in)
    if (!identical(ex_or_in_levels, EX_OR_IN_LEVELS2)
     && !identical(ex_or_in_levels, EX_OR_IN_LEVELS))
        stop("\"ex_or_in\" column has invalid levels")
    if (!("label.color" %in% extract_colnames))
        ans$label.color <- "blue"
    if (!("lty" %in% extract_colnames))
        ans$lty <- c("solid", "solid", "dashed", "solid")[ex_or_in]
    if (!("color" %in% extract_colnames))
        ans$color <- c("orange", "darkgrey", "grey", "black")[ex_or_in]
    if (!("width" %in% extract_colnames)
     && "UATXHcount" %in% extract_colnames) {
        min_UATXHcount <- min(ans$UATXHcount)
        if (min_UATXHcount < 0L) {
            warning("'UATXHcount' column contains negative values. Cannot use ",
                    "it to set the widths of the edges.")
        } else {
            max_UATXHcount <- max(ans$UATXHcount)
            if (max_UATXHcount <= 0L) {
                warning("'UATXHcount' column has no positive values. Cannot use ",
                        "it to set the widths of the edges.")
            } else {
                ans$width <- 20.0 * ans$UATXHcount / max(ans$UATXHcount)
            }
        }
    }
    ans
}

.make_igraph <- function(d)
{
    ## Prepare the 'vertices' argument to pass to graph.data.frame().
    nodes <- sgnodes(d)
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
    #g$layout <- layout.sgraph

    g
}

### 'sgedges0' must be a data.frame as returned by:
###     sgedges( , keep.dup.edges=TRUE)
.make_igraph_from_sgedges0 <- function(sgedges0, gene_id=NA,
                                       tx_id.as.edge.label=FALSE)
{
    if (!is.data.frame(sgedges0))
        stop("'sgedges0' must be a data.frame")
    if (!isTRUEorFALSE(tx_id.as.edge.label))
        stop("'tx_id.as.edge.label' must be TRUE or FALSE")
    d <- .precook_igraph_edges_from_sgedges(sgedges0)
    if (tx_id.as.edge.label)
        d$label <- d$tx_id
    .make_igraph(d)
}

### 'sgedges' must be a DataFrame as returned by:
###     sgedges( , keep.dup.edges=FALSE)
### or by:
###     sgedges2( )
.make_igraph_from_sgedges <- function(sgedges, gene_id=NA,
                                      tx_id.as.edge.label=FALSE)
{
    if (!is(sgedges, "DataFrame"))
        stop("'sgedges' must be a DataFrame")
    if (!isTRUEorFALSE(tx_id.as.edge.label))
        stop("'tx_id.as.edge.label' must be TRUE or FALSE")
    d <- .precook_igraph_edges_from_sgedges(sgedges)
    if (tx_id.as.edge.label)
        d$label <- sapply(d$tx_id, paste, collapse=",")
    d$tx_id <- NULL
    ## Turning 'd' into an ordinary data.frame. (Looks like 'as.data.frame()'
    ## on a DataFrame ignores the 'stringsAsFactors' arg so we use
    ## 'data.frame(as.list())' instead.)
    d <- data.frame(as.list(d), stringsAsFactors=FALSE)
    .make_igraph(d)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgraph() extractor
###
### Returns the splicing graph in an Ragraph object.
###

setGeneric("sgraph", signature="x",
    function(x, gene_id=NA, keep.dup.edges=FALSE,
             tx_id.as.edge.label=FALSE, as.igraph=FALSE)
        standardGeneric("sgraph")
)

setMethod("sgraph", "ANY",
    function(x, gene_id=NA, keep.dup.edges=FALSE,
             tx_id.as.edge.label=FALSE, as.igraph=FALSE)
    {
        sgedges <- sgedges(x, gene_id=gene_id, keep.dup.edges=keep.dup.edges)
        sgraph(sgedges, tx_id.as.edge.label=tx_id.as.edge.label,
                        as.igraph=as.igraph)
    }
)

setMethod("sgraph", "data.frame",
    function(x, gene_id=NA, keep.dup.edges=FALSE,
             tx_id.as.edge.label=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a data.frame")
        igraph <- .make_igraph_from_sgedges0(x,
                      tx_id.as.edge.label=tx_id.as.edge.label)
        sgraph(igraph, as.igraph=as.igraph)
    }
)

setMethod("sgraph", "DataFrame",
    function(x, gene_id=NA, keep.dup.edges=FALSE,
             tx_id.as.edge.label=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a DataFrame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a DataFrame")
        igraph <- .make_igraph_from_sgedges(x,
                      tx_id.as.edge.label=tx_id.as.edge.label)
        sgraph(igraph, as.igraph=as.igraph)
    }
)

setMethod("sgraph", "igraph",
    function(x, gene_id=NA, keep.dup.edges=FALSE,
             tx_id.as.edge.label=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an igraph object")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is an igraph object")
        if (!identical(tx_id.as.edge.label, FALSE))
            stop("the 'tx_id.as.edge.label' arg is not supported ",
                 "when 'x' is an igraph object")
        if (!isTRUEorFALSE(as.igraph))
            stop("'as.igraph' must be TRUE or FALSE")
        if (as.igraph) {
            ## Need to load the igraph package so the user can display, plot,
            ## and manipulate the returned object.
            library(igraph)
            return(x)  # no-op
        }
        make_Ragraph_from_igraph(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sgraph2() extractor
###
### Same as sgraph() except that uninformative nodes (i.e. SSids) are removed.
###

sgraph2 <- function(x, gene_id=NA, tx_id.as.edge.label=FALSE, as.igraph=FALSE)
{
    sgraph(sgedges2(x, gene_id=gene_id),
           tx_id.as.edge.label=tx_id.as.edge.label, as.igraph=as.igraph)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "plot" method.
###

setMethod("plot", c("SplicingGraphs", "ANY"),
    function(x, y, gene_id=NA)
    {
        if (missing(gene_id)) {
            if (missing(y)) {
                gene_id <- NA
            } else {
                gene_id <- y
            }
        } else {
            if (!missing(y))
                warning("'y' is ignored when plotting a SplicingGraphs ",
                        "object and 'gene_id' is supplied")
        }
        if (!isSingleStringOrNA(gene_id))
            stop("the supplied gene id must be a single string (or NA)")
        x_names <- names(x)
        if (!is.null(x_names) && is.na(gene_id))
            stop("You need to specify a gene id when 'x' has names ",
                 "e.g. 'plot(sg, \"some gene id\")'. Get all valid ",
                 "gene ids with 'unique(names(sg))'.")
        plot(sgraph(x, gene_id=gene_id))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### slideshow()
###

slideshow <- function(x)
{
    if (!is(x, "SplicingGraphs"))
        stop("'x' must be a SplicingGraphs object")
    for (gene_id in unique(names(x))) {
        ntx <- sum(names(x) == gene_id)
        cat("Plotting gene ", gene_id, " (", ntx, " transcripts). ", sep="")
        plot(x, gene_id)
        cat("Press <Enter> for next...")
        readLines(n=1)
    }
}

