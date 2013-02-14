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

### Should return a CompressedIntegerList.
setMethod("Spath", "GRangesList",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        ans <- mcols(x)[ , "Spath"]
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
### UATXHcount() accessor
###

setGeneric("UATXHcount", signature="x",
    function(x, gene_id=NA) standardGeneric("UATXHcount")
)

### Should return an integer vector or a NULL.
setMethod("UATXHcount", "GRangesList",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        ans <- mcols(x)[["UATXHcount"]]
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
### .hits() accessor (not exported)
###

setGeneric(".hits", signature="x",
    function(x, gene_id=NA) standardGeneric(".hits")
)

### Should return a CompressedCharacterList or a NULL.
setMethod(".hits", "GRangesList",
    function(x, gene_id=NA)
    {
        if (!isSingleStringOrNA(gene_id))
            stop("'gene_id' must be a single string (or NA)")
        if (length(x) == 0L)
            stop("'x' must be of length >= 1")
        x_names <- names(x)
        if (is.null(x_names)) {
            if (!is.na(gene_id))
                stop("the 'gene_id' arg is not supported ",
                     "when 'x' is unnamed (in which case all its elements ",
                     "(i.e. transcripts) are considered to belong to the ",
                     "same gene)")
            ans <- mcols(unlist(x, use.names=FALSE))[["hits"]]
            return(ans) 
        }
        if (is.na(gene_id))
            stop("'gene_id' must be supplied when 'x' has names")
        x <- x[x_names == gene_id]
        if (length(x) == 0L)
            stop("invalid 'gene_id'")
        ans <- mcols(unlist(x, use.names=FALSE))[["hits"]]
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
### Returns a 4-col (or 5-col if 'UATXHcount' is supplied) data.frame representing
### the splicing graph.
.make_Sgdf0_from_Spath <- function(spath, UATXHcount=NULL)
{
    if (!is.null(UATXHcount)) {
        if (!is.integer(UATXHcount))
            stop("'UATXHcount' must be an integer vector or NULL")
        if (length(UATXHcount) != length(spath))
            stop("when not NULL, 'UATXHcount' must have ",
                 "the same length as 'spath'")
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
    if (!is.null(UATXHcount))
        sgdf0$UATXHcount <- rep.int(UATXHcount, nedges_per_tx)
    sgdf0
}

### Collapse the duplicated edges in 'sgdf0' into a DataFrame.
### We use a DataFrame instead of a data.frame because we want to store
### the tx_id col in a CompressedFactorList (even though this container
### doesn't formally exist and a CompressedIntegerList is actually what's
### being used).
.make_Sgdf_from_Sgdf0 <- function(sgdf0, ex_hits=NULL, in_hits=NULL)
{
    from <- sgdf0[ , "from"]
    to <- sgdf0[ , "to"]
    ex_or_in <- sgdf0[ , "ex_or_in"]
    tx_id <- sgdf0[ , "tx_id"]
    edges <- paste(from, to, sep="~")
    sm <- match(edges, edges)
    if (!all(ex_or_in == ex_or_in[sm]))
        stop("invalid splicing graph")
    is_not_dup <- sm == seq_along(sm)
    sgdf <- DataFrame(sgdf0[is_not_dup, , drop=FALSE])
    sgdf$tx_id <- splitAsList(tx_id, sm)
    UATXHcount <- sgdf$UATXHcount
    if (!is.null(UATXHcount))
        sgdf$UATXHcount <- sum(splitAsList(sgdf0$UATXHcount, sm))
    if (is.null(ex_hits) && is.null(in_hits))
        return(sgdf)
    hits <- relist(character(0), PartitioningByEnd(NG=length(sm)))
    if (!is.null(ex_hits)) {
        if (!is(ex_hits, "CharacterList"))
            stop("'ex_hits' must be a CharacterList object")
        ex_idx <- which(ex_or_in == "ex")
        if (length(ex_idx) != length(ex_hits))
            stop("'ex_hits' is incompatible with 'sgdf0'")
        hits[ex_idx] <- ex_hits
    }
    if (!is.null(in_hits)) {
        if (!is(in_hits, "CharacterList"))
            stop("'in_hits' must be a CharacterList object")
        in_idx <- which(ex_or_in == "in")
        if (length(in_idx) != length(in_hits))
            stop("'in_hits' is incompatible with 'sgdf0'")
        hits[in_idx] <- in_hits
    }
    ## TODO: This is quite inefficient. Improve it.
    for (i in which(!is_not_dup))
        hits[[sm[i]]] <- unique(hits[[sm[i]]], hits[[i]])
    sgdf$hits <- hits[is_not_dup]
    sgdf$nhits <- elementLengths(sgdf$hits)
    sgdf
}

setGeneric("Sgdf", signature="x",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
        standardGeneric("Sgdf")
)

setMethod("Sgdf", "ANY",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
    {
        spath <- Spath(x, gene_id=gene_id)
        if (is.null(UATXHcount))
            UATXHcount <- UATXHcount(x, gene_id=gene_id)
        if (is.null(inbytx))
            return(Sgdf(spath, UATXHcount=UATXHcount,
                               keep.dup.edges=keep.dup.edges))
        if (!is(inbytx, "GRangesList"))
            stop("'inbytx' must be NULL or a GRangesList object")
        if (!is(x, "GRangesList"))
            stop("'x' must be a GRangesList object ",
                 "when 'inbytx' is a GRangesList object")
        if (length(inbytx) != length(x))
            stop("'inbytx' must have the same length as 'x'")
        if (!identical(elementLengths(inbytx) + 1L, elementLengths(x)))
            stop("the shape of 'inbytx' is not compatible ",
                 "with the shape of 'x'")
        if (!identical(keep.dup.edges, FALSE))
            stop("'keep.dup.edges' must be FALSE when 'inbytx' is supplied")
        sgdf0 <- Sgdf(spath, UATXHcount=UATXHcount, keep.dup.edges=TRUE)
        ex_or_in <- sgdf0[ , "ex_or_in"]
        ex_hits <- .hits(x, gene_id=gene_id)
        if (is.null(ex_hits))
            stop("'x' must have a \"hits\" inner metadata column ",
                 "when 'inbytx' is a GRangesList object. May be ",
                 "you forgot to pass it thru assignSubfeatureHits()?")
        in_hits <- .hits(inbytx, gene_id=gene_id)
        if (is.null(in_hits))
            stop("'inbytx' has no \"hits\" inner metadata column. May be ",
                 "you forgot to pass it thru assignSubfeatureHits()?")
        .make_Sgdf_from_Sgdf0(sgdf0, ex_hits=ex_hits, in_hits=in_hits)
    }
)

setMethod("Sgdf", "IntegerList",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an IntegerList")
        if (!is.null(inbytx))
            stop("the 'inbytx' arg is not supported ",
                 "when 'x' is an IntegerList")
        sgdf0 <- .make_Sgdf0_from_Spath(x, UATXHcount=UATXHcount)
        Sgdf(sgdf0, keep.dup.edges=keep.dup.edges)
    }
)

setMethod("Sgdf", "data.frame",
    function(x, gene_id=NA, UATXHcount=NULL, inbytx=NULL, keep.dup.edges=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!is.null(UATXHcount))
            stop("the 'UATXHcount' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!is.null(inbytx))
            stop("the 'inbytx' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!isTRUEorFALSE(keep.dup.edges))
            stop("'keep.dup.edges' must be TRUE or FALSE")
        if (keep.dup.edges)
            return(x)  # no-op
        .make_Sgdf_from_Sgdf0(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_igraph_from_Sgdf()
###

### 'sgdf' must be a data.frame as returned by:
###     Sgdf( , keep.dup.edges=TRUE)
### or a DataFrame as returned by:
###     Sgdf( , keep.dup.edges=FALSE)
### Valid extra cols are: "label", "label.color", "lty", "color", "width"
### and "UATXHcount". They are used to set graphical parameters on the edges.
.precook_igraph_edges_from_Sgdf <- function(sgdf)
{
    required_colnames <- c("from", "to", "ex_or_in", "tx_id")
    extra_colnames <- c("label", "label.color", "lty", "color",
                        "width", "UATXHcount")
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
    #g$layout <- layout.Sgraph

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sgraph() accessor
###
### Returns the splicing graph in an igraph object.
###

setGeneric("Sgraph", signature="x",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
        standardGeneric("Sgraph")
)

setMethod("Sgraph", "ANY",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        sgdf <- Sgdf(x, gene_id=gene_id, keep.dup.edges=keep.dup.edges)
        Sgraph(sgdf, as.igraph=as.igraph)
    }
)

setMethod("Sgraph", "data.frame",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a data.frame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a data.frame")
        igraph <- .make_igraph_from_Sgdf0(x)
        Sgraph(igraph, as.igraph=as.igraph)
    }
)

setMethod("Sgraph", "DataFrame",
    function(x, gene_id=NA, keep.dup.edges=FALSE, as.igraph=FALSE)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is a DataFrame")
        if (!identical(keep.dup.edges, FALSE))
            stop("the 'keep.dup.edges' arg is not supported ",
                 "when 'x' is a DataFrame")
        igraph <- .make_igraph_from_Sgdf(x)
        Sgraph(igraph, as.igraph=as.igraph)
    }
)

setMethod("Sgraph", "igraph",
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

Sgraph2 <- function(x, gene_id=NA, as.igraph=FALSE)
{
    Sgraph(Sgdf2(x, gene_id=gene_id), as.igraph=as.igraph)
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

