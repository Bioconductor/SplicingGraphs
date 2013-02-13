### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A layout function for igraph objects.
###

### Experimental. Not ready yet!
layout.Sgraph <- function(graph)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_Ragraph_from_igraph()
###

### Safer than using sub().
.safeTranslateEdgeNames <- function(edge_names, old.sep="|", new.sep="~")
{
    edge_names <- strsplit(edge_names, old.sep, fixed=TRUE)
    stopifnot(all(elementLengths(edge_names) == 2L))  # sanity check
    edge_names <- unlist(edge_names, use.names=FALSE)
    paste0(edge_names[c(TRUE, FALSE)], new.sep, edge_names[c(FALSE, TRUE)])
}

### Mappings from igraph attribute names to Rgraphviz attribute names.
.IGRAPH_2_RGRAPHVIZ_EDGE_ATTRNAMES <- c(
    width="lwd",
    label.color="fontcolor"
    #lty="style"
)

.IGRAPH_2_RGRAPHVIZ_NODE_ATTRNAMES <- c(
    color="fillcolor",
    .IGRAPH_2_RGRAPHVIZ_EDGE_ATTRNAMES
)

.safeTranslateAttrNames <- function(attr_names, old2new)
{
    m <- match(attr_names, names(old2new), nomatch=0L)
    idx <- which(m != 0L)
    new_names <- old2new[m]
    do_it <- !(new_names %in% attr_names)
    attr_names[idx[do_it]] <- new_names[do_it]
    attr_names
}

.make_Ragraph_nodeAttrs_from_graphNEL <- function(graph_nel)
{
    attr_names <- names(nodeDataDefaults(graph_nel))
    node_data <- nodeData(graph_nel)
    node_names <- names(node_data)
    node_attrs <- lapply(attr_names,
                         function(attr_name)
                             unlist(lapply(node_data, `[[`, attr_name),
                                    recursive=FALSE))
    names(node_attrs) <- .safeTranslateAttrNames(attr_names,
                                 .IGRAPH_2_RGRAPHVIZ_NODE_ATTRNAMES)
    is_null <- sapply(node_attrs, is.null)
    node_attrs <- node_attrs[!is_null]

    node_fontsize <- rep.int("12", length(node_names))
    names(node_fontsize) <- node_names
    node_attrs$fontsize <- node_fontsize

    node_attrs
}

.make_Ragraph_edgeAttrs_from_graphNEL <- function(graph_nel)
{
    attr_names <- names(edgeDataDefaults(graph_nel))
    edge_data <- edgeData(graph_nel)
    edge_data_names <- names(edge_data)
    if (length(attr_names) == 0L || length(edge_data_names) == 0L)
        return(list())
    if (anyDuplicated(edge_data_names))
        warning("graph object has more than 1 edge between the same ",
                "2 nodes. Because setting plotting attributes for those ",
                "edges (e.g. label, color, line width, line style, etc...) ",
                "is not fully supported yet, they will end up with the same ",
                "attributes. Which is unlikely to be what you want.")
    edge_data_names <- .safeTranslateEdgeNames(edge_data_names)

    ## edgeNames() is broken on graphNEL objects that have more than 1 edge
    ## between the same 2 nodes, hence the use of unique() in the sanity check.
    edge_names <- edgeNames(graph_nel)
    stopifnot(identical(unique(edge_data_names), edge_names))  # sanity check

    names(edge_data) <- edge_data_names
    edge_attrs <- lapply(attr_names,
                         function(attr_name)
                             unlist(lapply(edge_data, `[[`, attr_name),
                                    recursive=FALSE))
    names(edge_attrs) <- .safeTranslateAttrNames(attr_names,
                                 .IGRAPH_2_RGRAPHVIZ_EDGE_ATTRNAMES)
    is_null <- sapply(edge_attrs, is.null)
    edge_attrs <- edge_attrs[!is_null]

    edge_fontsize <- rep.int("10", length(edge_data_names))
    names(edge_fontsize) <- edge_data_names
    edge_attrs$fontsize <- edge_fontsize

    edge_attrs
}

.make_Ragraph_from_graphNEL <- function(graph_nel, gene_id=NA)
{
    node_attrs <- .make_Ragraph_nodeAttrs_from_graphNEL(graph_nel)
    edge_attrs <- .make_Ragraph_edgeAttrs_from_graphNEL(graph_nel)
    agopen(graph_nel, name=gene_id, nodeAttrs=node_attrs, edgeAttrs=edge_attrs)
}

make_Ragraph_from_igraph <- function(igraph, gene_id=NA)
{
    .make_Ragraph_from_graphNEL(igraph.to.graphNEL(igraph), gene_id=gene_id)
}

