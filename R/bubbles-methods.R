### =========================================================================
### "bubbles" methods
### -------------------------------------------------------------------------


### Returns a DataFrame with 1 row per bubble and the following cols:
### source <character>, sink <character>, left_path <IntegerList>, and
### right_path <IntegerList>.
### Must find at least 1 bubble if 'left_txpath' and 'right_txpath' are not
### identical.
.find_bubbles_in_tx_pair <- function(left_txpath, right_txpath)
{
    ans_source <- ans_sink <- character(0)
    ans_left_path <- ans_right_path <- IntegerList()
    left_len <- length(left_txpath)
    right_len <- length(right_txpath)
    i <- j <- 1L
    in_bubble <- FALSE
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
            if (in_bubble) {
                ## End of the current bubble.
                in_bubble <- FALSE
                if (curr_left_SSid == .Machine$integer.max) {
                    sink <- "L"
                } else {
                    sink <- as.character(curr_left_SSid)
                }
                ans_sink <- c(ans_sink, sink)
                ans_left_path <- c(ans_left_path,
                                   IntegerList(bubble_left_path))
                ans_right_path <- c(ans_right_path,
                                    IntegerList(bubble_right_path))
            }
            i <- i + 1L
            j <- j + 1L
            next
        }
        if (!in_bubble) {
            ## Start a new bubble.
            in_bubble <- TRUE
            if (i == 1L) {
                bubble_source <- "R"
            } else {
                bubble_source <- as.character(left_txpath[i - 1L])
            }
            ans_source <- c(ans_source, bubble_source)
            bubble_left_path <- bubble_right_path <- integer(0)
        }
        if (curr_left_SSid < curr_right_SSid) {
            bubble_left_path <- c(bubble_left_path, curr_left_SSid)
            i <- i + 1L
        } else {
            bubble_right_path <- c(bubble_right_path, curr_right_SSid)
            j <- j + 1L
        }
    }
    DataFrame(source=ans_source,
              sink=ans_sink,
              left_path=ans_left_path,
              right_path=ans_right_path)
}

.extract_bubbles_from_txpaths <- function(txpaths)
{
    ntx <- length(txpaths)
    if (ntx <= 1L) {
        ## No bubbles.
        bubbles <- .find_bubbles_in_tx_pair(integer(0), integer(0))
        bubbles <- cbind(bubbles, DataFrame(left_tx_id=character(0),
                                            right_tx_id=character(0)))
        return(bubbles)
    }
    tx_id <- names(txpaths)
    if (is.null(tx_id))
        tx_id <- seq_along(txpaths)

    npairs <- (ntx * (ntx - 1L)) %/% 2L
    all_bubbles <- vector(mode="list", length=npairs)
    z <- 1L
    for (i in 1:(ntx-1L)) {
        left_txpath <- txpaths[[i]]
        left_tx_id <- tx_id[i]
        for (j in (i+1L):ntx) {
            right_txpath <- txpaths[[j]]
            right_tx_id <- tx_id[j]
            bubbles <- .find_bubbles_in_tx_pair(left_txpath, right_txpath)
            bubbles <- cbind(bubbles, DataFrame(left_tx_id=left_tx_id,
                                                right_tx_id=right_tx_id))
            all_bubbles[[z]] <- bubbles
            z <- z + 1L
        }
    }
    ans <- do.call(rbind, all_bubbles)

    ## Remove duplicate bubbles.
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
        .extract_bubbles_from_txpaths(x)
    }
)

### TODO: Add "bubbles" methods for data.frame and DataFrame objects.

