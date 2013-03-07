### =========================================================================
### "bubbles" methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_ASevents2_from_txpaths()
###

### Returns a DataFrame with 1 row per AS event and the following cols:
### source <character>, sink <character>, left_path <IntegerList>, and
### right_path <IntegerList>.
### Must find at least 1 AS event if 'left_txpath' and 'right_txpath' are
### not identical.
.extract_ASevents2_from_tx_pair <- function(left_txpath, right_txpath)
{
    ans_source <- ans_sink <- character(0)
    ans_left_path <- ans_right_path <- IntegerList()
    left_len <- length(left_txpath)
    right_len <- length(right_txpath)
    i <- j <- 1L
    in_ASevent <- FALSE
    while (i <= left_len + 1L || j <= right_len + 1L) {
        if (i <= left_len) {
            curr_left_SSid <- left_txpath[i]
        } else {
            curr_left_SSid <- .Machine$integer.max
        }
        if (j <= right_len) {
            curr_right_SSid <- right_txpath[j]
        } else {
            curr_right_SSid <- .Machine$integer.max
        }
        if (curr_left_SSid == curr_right_SSid) {
            if (in_ASevent) {
                ## End of the current AS event.
                in_ASevent <- FALSE
                if (curr_left_SSid == .Machine$integer.max) {
                    sink <- "L"
                } else {
                    sink <- as.character(curr_left_SSid)
                }
                ans_sink <- c(ans_sink, sink)
                ans_left_path <- c(ans_left_path,
                                   IntegerList(ASevent_left_path))
                ans_right_path <- c(ans_right_path,
                                    IntegerList(ASevent_right_path))
            }
            i <- i + 1L
            j <- j + 1L
            next
        }
        if (!in_ASevent) {
            ## Start a new AS event.
            in_ASevent <- TRUE
            if (i == 1L) {
                ASevent_source <- "R"
            } else {
                ASevent_source <- as.character(left_txpath[i - 1L])
            }
            ans_source <- c(ans_source, ASevent_source)
            ASevent_left_path <- ASevent_right_path <- integer(0)
        }
        if (curr_left_SSid < curr_right_SSid) {
            ASevent_left_path <- c(ASevent_left_path, curr_left_SSid)
            i <- i + 1L
        } else {
            ASevent_right_path <- c(ASevent_right_path, curr_right_SSid)
            j <- j + 1L
        }
    }
    DataFrame(source=ans_source,
              sink=ans_sink,
              left_path=ans_left_path,
              right_path=ans_right_path)
}

.extract_ASevents2_from_txpaths <- function(txpaths)
{
    ntx <- length(txpaths)
    if (ntx <= 1L) {
        ## No AS events.
        ASevents2 <- .extract_ASevents2_from_tx_pair(integer(0), integer(0))
        ans <- cbind(ASevents2, DataFrame(left_tx_id=character(0),
                                          right_tx_id=character(0)))
        return(ans)
    }
    tx_id <- names(txpaths)
    if (is.null(tx_id))
        tx_id <- seq_along(txpaths)

    npairs <- (ntx * (ntx - 1L)) %/% 2L
    all_ASevents2 <- vector(mode="list", length=npairs)
    z <- 1L
    for (i in 1:(ntx-1L)) {
        left_txpath <- txpaths[[i]]
        left_tx_id <- tx_id[i]
        for (j in (i+1L):ntx) {
            right_txpath <- txpaths[[j]]
            right_tx_id <- tx_id[j]
            ASevents2 <- .extract_ASevents2_from_tx_pair(left_txpath,
                                                         right_txpath)
            ASevents2 <- cbind(ASevents2, DataFrame(left_tx_id=left_tx_id,
                                                    right_tx_id=right_tx_id))
            all_ASevents2[[z]] <- ASevents2
            z <- z + 1L
        }
    }
    ans <- do.call(rbind, all_ASevents2)

    ## Remove duplicate AS events.
    keys <- c("source", "sink", "left_path", "right_path")
    key1 <- ans[ , "source"]
    key2 <- ans[ , "sink"]
    key3 <- unlist(lapply(ans[ , "left_path"], paste, collapse=","),
                   use.names=FALSE)
    key4 <- unlist(lapply(ans[ , "right_path"], paste, collapse=","),
                   use.names=FALSE)
    key <- paste(key1, key2, key3, key4, sep="|")
    left_tx_id <- ans[ , "left_tx_id"]
    right_tx_id <- ans[ , "right_tx_id"]
    sm <- match(key, key)
    is_not_dup <- sm == seq_along(sm)
    ans <- ans[is_not_dup, keys, drop=FALSE]
    rownames(ans) <- NULL
    ans$left_tx_id <- unique(splitAsList(left_tx_id, sm))
    ans$right_tx_id <- unique(splitAsList(right_tx_id, sm))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_bubbles_from_txpathmat()
###

### Returns a DataFrame with 1 row per bubble and the following cols:
### source <character>, sink <character>, d <integer>.
.extract_bubbles_from_txpathmat <- function(txpathmat)
{
    stop("not ready yet, sorry!")
    ncol <- ncol(txpathmat)
    for (i in 1:(ncol-2L)) {
        for (j in (i+2L):ncol) {
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bubbles()
###

setGeneric("bubbles", signature="x",
    function(x, gene_id=NA) standardGeneric("bubbles")
)

setMethod("bubbles", "ANY",
    function(x, gene_id=NA)
    {
        txpaths <- txpaths(x, gene_id=gene_id)
        bubbles(txpaths)
    }
)

setMethod("bubbles", "IntegerList",
    function(x, gene_id=NA)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an IntegerList")
        #.extract_ASevents2_from_txpaths(x)
        txpathmat <- make_matrix_from_txpaths(x)
        .extract_bubbles_from_txpathmat(txpathmat)
    }
)

### TODO: Add "bubbles" methods for data.frame and DataFrame objects.

