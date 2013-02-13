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
    width="weight",
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
    edge_names <- edgeNames(graph_nel)
    if (length(attr_names) == 0L || length(edge_names) == 0L)
        return(list())
    if (anyDuplicated(edge_names))
        warning("graphNEL object has more than 1 edge between the same ",
                "2 nodes. Plotting attributes for those edges (e.g. label, ",
                "color, line width, line style, etc...) are likely to be ",
                "wrong.")
    edge_data <- edgeData(graph_nel)
    stopifnot(identical(.safeTranslateEdgeNames(names(edge_data)),
                        edge_names))  # sanity check
    names(edge_data) <- edge_names
    edge_attrs <- lapply(attr_names,
                         function(attr_name)
                             unlist(lapply(edge_data, `[[`, attr_name),
                                    recursive=FALSE))
    names(edge_attrs) <- .safeTranslateAttrNames(attr_names,
                                 .IGRAPH_2_RGRAPHVIZ_EDGE_ATTRNAMES)
    is_null <- sapply(edge_attrs, is.null)
    edge_attrs <- edge_attrs[!is_null]

    edge_fontsize <- rep.int("10", length(edge_names))
    names(edge_fontsize) <- edge_names
    edge_attrs$fontsize <- edge_fontsize

    edge_attrs
}

make_Ragraph_from_graphNEL <- function(graph_nel, gene_id=NA)
{
    node_attrs <- .make_Ragraph_nodeAttrs_from_graphNEL(graph_nel)
    edge_attrs <- .make_Ragraph_edgeAttrs_from_graphNEL(graph_nel)
    agopen(graph_nel, name=gene_id, nodeAttrs=node_attrs, edgeAttrs=edge_attrs)
}

