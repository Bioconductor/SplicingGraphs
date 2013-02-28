### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A layout function for a splicing graph represented as an igraph object.
###

### Experimental. Not ready yet!
layout.sgraph <- function(graph)
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
### .igraph.to.graphNEL2() - Replacement for igraph::igraph.to.graphNEL()
###
### igraph::igraph.to.graphNEL() works on a graph with more than 1 edge
### between the same 2 nodes: surprisingly it returns a graphNEL object that
### is plottable, seems to have the correct number of edges, and passes the
### validGraph() test. Nonetheless it's broken! In particular, the edge
### attributes (aka the "edge data" in graphNEL jargon, found in
### x@edgeData@data for graphNEL object 'x') propagated from the original
### igraph object are corrupted. The .igraph.to.graphNEL2() function below
### fixes that, and only that (i.e. it repairs x@edgeData@data), because this
### is all what we care about.
###

.find_old2new_mapping <- function(old, new)
{
    N <- length(old)
    if (length(new) != N)
        stop("'old' and 'new' must have the same length")
    oo1 <- order(old)
    oo2 <- order(new)
    ans <- integer(N)
    ans[oo2] <- oo1
    stopifnot(all(old[ans] == new))
    ans
}

.repair_edge_data <- function(edge_data, graph)
{
    old_edges <- get.data.frame(graph)
    stopifnot(identical(colnames(old_edges)[1:2], c("from", "to")))
    old_edge_names <- paste0(old_edges[ , "from"], "|", old_edges[ , "to"])
    new_edge_names <- names(edge_data)
    old2new <- .find_old2new_mapping(old_edge_names, new_edge_names)
    old_edges <- old_edges[old2new, -(1:2), drop=FALSE]
    edge_data <- lapply(seq_along(edge_data),
                            function(i) {
                                old_data <- as.list(old_edges[i, ])
                                new_data <- edge_data[[i]]
                                new_data[names(old_data)] <- old_data
                                new_data
                            })
    names(edge_data) <- new_edge_names
    edge_data
}

.igraph.to.graphNEL2 <- function(graph)
{
    ans <- igraph.to.graphNEL(graph)
    ## We use direct slot access instead of the edgeData() accessor (which
    ## is broken on a graph with more than 1 edge between the same 2 nodes).
    ans@edgeData@data <- .repair_edge_data(ans@edgeData@data, graph)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_Ragraph_from_igraph()
###

### Mappings from igraph attribute names to Ragraph attribute names.
.igraph2Ragraph_EDGE_ATTRNAMES <- c(
    width="lwd",
    label.color="fontcolor"
    #lty="style"
)

.igraph2Ragraph_NODE_ATTRNAMES <- c(
    color="fillcolor",
    .igraph2Ragraph_EDGE_ATTRNAMES
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

### Safer than using sub().
.safeTranslateEdgeNames <- function(edge_names, old.sep="|", new.sep="~")
{
    edge_names <- strsplit(edge_names, old.sep, fixed=TRUE)
    stopifnot(all(elementLengths(edge_names) == 2L))  # sanity check
    edge_names <- unlist(edge_names, use.names=FALSE)
    paste0(edge_names[c(TRUE, FALSE)], new.sep, edge_names[c(FALSE, TRUE)])
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
                                 .igraph2Ragraph_NODE_ATTRNAMES)
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
    ## We use direct slot access instead of the edgeData() accessor (which
    ## is broken on a graph with more than 1 edge between the same 2 nodes).
    edge_data <- graph_nel@edgeData@data
    edge_data_names <- names(edge_data)
    if (length(attr_names) == 0L || length(edge_data_names) == 0L)
        return(list())
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
                                 .igraph2Ragraph_EDGE_ATTRNAMES)
    is_null <- sapply(edge_attrs, is.null)
    edge_attrs <- edge_attrs[!is_null]

    edge_fontsize <- rep.int("10", length(edge_data_names))
    names(edge_fontsize) <- edge_data_names
    edge_attrs$fontsize <- edge_fontsize

    edge_attrs
}

### Replacement for Rgraphviz::buildEdgeList().
### Rgraphviz::buildEdgeList() is broken when some of the components in the
### 'edgeAttrs' list have duplicated edge names on them (which only happens
### when the graph has more than 1 edge between the same 2 nodes).
### Note that .buildEdgeList2() assumes that 'edgeAttrs' was extracted from
### 'graph' with .make_Ragraph_edgeAttrs_from_graphNEL(), and that 'graph'
### still contains the original edge data. Because of this, .buildEdgeList2()
### is not a serious candidate for replacing Rgraphviz::buildEdgeList().

.get_pEdge_fixed_attrs <- function(pedge_attrs, edge_names, edge_pos, edgeAttrs)
{
    edge_name <- edge_names[edge_pos]
    ## Check and fix 'attrs' slot.
    attrs0 <- sapply(edgeAttrs,
                     function(edge_attr) {
                         idx <- which(names(edge_attr) == edge_name)
                         if (length(idx) == 0L)
                             return(NA)
                         if (length(idx) == 1L)
                             return(edge_attr[[idx]])
                         if (!identical(names(edge_attr), edge_names))
                             stop("cannot fix pEdge nb ", edge_pos)
                         edge_attr[[edge_pos]]
                     })
    if (!is.character(attrs0))
        storage.mode(attrs0) <- attrs0
    attrs0 <- attrs0[!is.na(attrs0)]
    attrs0_names <- names(attrs0)
    m <- match(attrs0_names, names(pedge_attrs))
    stopifnot(!any(is.na(m)))
    attrs1 <- unlist(pedge_attrs[attrs0_names], recursive=FALSE)
    as.list(attrs0[attrs1 != attrs0])
}

.buildEdgeList2 <- function(graph, recipEdges=c("combined", "distinct"),
                            edgeAttrs=list(), subGList=list(),
                            defAttrs=list())
{
    ans <- buildEdgeList(graph, recipEdges=recipEdges,
                         edgeAttrs=edgeAttrs, subGList=subGList,
                         defAttrs=defAttrs)
    ans_names <- names(ans)  # the edge names
    ## Checking that 'graph' still contains the original edge data.
    ## We use direct slot access instead of the edgeData() accessor (which
    ## is broken on a graph with more than 1 edge between the same 2 nodes).
    stopifnot(identical(ans_names,
              .safeTranslateEdgeNames(names(graph@edgeData@data))))
    ans <- lapply(seq_along(ans),
                  function(i) {
                      pedge <- ans[[i]]
                      ## Check 'from' and 'to' slots.
                      edge_name <- paste0(from(pedge), "~", to(pedge))
                      stopifnot(identical(edge_name, ans_names[i]))
                      ## Get the fixed attribs (as a named list).
                      fixed_attrs <- .get_pEdge_fixed_attrs(pedge@attrs,
                                                            ans_names, i,
                                                            edgeAttrs)
                      ## Set them.
                      pedge@attrs[names(fixed_attrs)] <- fixed_attrs
                      pedge
                  })
    names(ans) <- ans_names
    ans
}

.make_Ragraph_from_graphNEL <- function(graph_nel, gene_id=NA)
{
    node_attrs <- .make_Ragraph_nodeAttrs_from_graphNEL(graph_nel)
    edge_attrs <- .make_Ragraph_edgeAttrs_from_graphNEL(graph_nel)
    ## Passing 'edge_attrs' directly to agopen() (thru the 'edgeAttrs' arg) is
    ## broken when there is more than 1 edge between the same 2 nodes, so we
    ## need to build and pass the lists of pNode and pEdge objects instead.
    #agopen(graph_nel, name=gene_id, nodeAttrs=node_attrs, edgeAttrs=edge_attrs)
    nodes <- buildNodeList(graph_nel, nodeAttrs=node_attrs)
    edges <- .buildEdgeList2(graph_nel, edgeAttrs=edge_attrs)
    agopen(name=gene_id, nodes=nodes, edges=edges,
                         edgeMode=edgemode(graph_nel))
}

make_Ragraph_from_igraph <- function(igraph, gene_id=NA)
{
    .make_Ragraph_from_graphNEL(.igraph.to.graphNEL2(igraph), gene_id=gene_id)
}

