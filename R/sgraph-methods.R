### =========================================================================
### sgraph (and related) methods
### -------------------------------------------------------------------------


setOldClass("igraph")

.EDGE_WEIGHTS <- c(1, 0.2, 0.1, 0.4)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_igraph_from_sgdf()
###

### 'sgdf' must be a data.frame as returned by:
###     sgdf( , keep.dup.edges=TRUE)
### or a DataFrame as returned by:
###     sgdf( , keep.dup.edges=FALSE)
### Valid extra cols are: "label", "label.color", "lty", "color", "width"
### and "UATXHcount". They are used to set graphical parameters on the edges.
.precook_igraph_edges_from_sgdf <- function(sgdf)
{
    required_colnames <- c("from", "to", "ex_or_in", "tx_id")
    extra_colnames <- c("label", "label.color", "lty", "color",
                        "width", "UATXHcount")
    extract_colnames <- c(required_colnames,
                          intersect(extra_colnames, colnames(sgdf)))
    ans <- sgdf[ , extract_colnames, drop=FALSE]
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
        ans$color <- c("green3", "darkgrey", "grey", "black")[ex_or_in]
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
    #g$layout <- layout.sgraph

    g
}

### 'sgdf0' must be a data.frame as returned by:
###     sgdf( , keep.dup.edges=TRUE)
.make_igraph_from_sgdf0 <- function(sgdf0, gene_id=NA)
{
    if (!is.data.frame(sgdf0))
        stop("'sgdf0' must be a data.frame")
    d <- .precook_igraph_edges_from_sgdf(sgdf0)
    if (!("label" %in% colnames(d)))
        d$label <- d$tx_id
    .make_igraph(d)
}

### 'sgdf' must be a DataFrame as returned by:
###     sgdf( , keep.dup.edges=FALSE)
### or by:
###     sgdf2( )
.make_igraph_from_sgdf <- function(sgdf, gene_id=NA)
{
    if (!is(sgdf, "DataFrame"))
        stop("'sgdf' must be a DataFrame")
    d <- .precook_igraph_edges_from_sgdf(sgdf)
    if (!("label" %in% colnames(d)))
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
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
        standardGeneric("sgraph")
)

setMethod("sgraph", "ANY",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        sgdf <- sgdf(x, gene_id=gene_id, keep.dup.edges=keep.dup.edges)
        sgraph(sgdf, as.igraph=as.igraph)
    }
)

setMethod("sgraph", "data.frame",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a data.frame")
        igraph <- .make_igraph_from_sgdf0(x)
        sgraph(igraph, as.igraph=as.igraph)
    }
)

setMethod("sgraph", "DataFrame",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a DataFrame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a DataFrame")
        igraph <- .make_igraph_from_sgdf(x)
        sgraph(igraph, as.igraph=as.igraph)
    }
)

setMethod("sgraph", "igraph",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an igraph object")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
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

sgraph2 <- function(x, gene_id=NA, as.igraph=FALSE)
{
    sgraph(sgdf2(x, gene_id=gene_id), as.igraph=as.igraph)
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

