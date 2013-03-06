### =========================================================================
### "bubbles" methods
### -------------------------------------------------------------------------


### Returns a DataFrame with 1 row per bubble and the following cols:
### source <character>, sink <character>, left_path <IntegerList>, and
### right_path <IntegerList>.
### Must find at least 1 bubble if 'left_SSids' and 'right_SSids' are not
### identical.
.find_bubbles_in_tx_pair <- function(left_SSids, right_SSids)
{
    ans_source <- ans_sink <- character(0)
    ans_left_path <- ans_right_path <- IntegerList()
    left_len <- length(left_SSids)
    right_len <- length(right_SSids)
    i <- j <- 1L
    prev_left <- prev_right <- "R"
    in_bubble <- FALSE
    while (i <= left_len + 1L || j <= right_len + 1L) {
        if (i <= left_len) {
            curr_left <- left_SSids[i]
        } else {
            curr_left <- .Machine$integer.max
        }
        if (j <= right_len) {
            curr_right <- right_SSids[j]
        } else {
            curr_right <- .Machine$integer.max
        }
        if (curr_left == curr_right) {
            if (in_bubble) {
                ## End of the current bubble.
                in_bubble <- FALSE
                if (curr_left == .Machine$integer.max) {
                    sink <- "L"
                } else {
                    sink <- as.character(curr_left)
                }
                ans_sink <- c(ans_sink, sink)
                ans_left_path <- c(ans_left_path,
                                   IntegerList(bubble_left_path))
                ans_right_path <- c(ans_right_path,
                                    IntegerList(bubble_right_path))
            }
            prev_left <- curr_left
            i <- i + 1L
            prev_right <- curr_right
            j <- j + 1L
            next
        }
        if (!in_bubble) {
            ## Start a new bubble.
            in_bubble <- TRUE
            ans_source <- c(ans_source, as.character(prev_left))
            bubble_left_path <- bubble_right_path <- integer(0)
        }
        if (curr_left < curr_right) {
            bubble_left_path <- c(bubble_left_path, curr_left)
            prev_left <- curr_left
            i <- i + 1L
        } else {
            bubble_right_path <- c(bubble_right_path, curr_right)
            prev_right <- curr_right
            j <- j + 1L
        }
    }
    DataFrame(source=ans_source,
              sink=ans_sink,
              left_path=ans_left_path,
              right_path=ans_right_path)
}

.extract_bubbles_from_spath <- function(spath)
{
    ntx <- length(spath)
    if (ntx <= 1L) {
        ## No bubbles.
        bubbles <- .find_bubbles_in_tx_pair(integer(0), integer(0))
        bubbles <- cbind(bubbles, DataFrame(left_tx_id=character(0),
                                            right_tx_id=character(0)))
        return(bubbles)
    }
    tx_id <- names(spath)
    if (is.null(tx_id))
        tx_id <- seq_along(spath)

    npairs <- (ntx * (ntx - 1L)) %/% 2L
    all_bubbles <- vector(mode="list", length=npairs)
    z <- 1L
    for (i in 1:(ntx-1L)) {
        left_SSids <- spath[[i]]
        left_tx_id <- tx_id[i]
        for (j in (i+1L):ntx) {
            right_SSids <- spath[[j]]
            right_tx_id <- tx_id[j]
            bubbles <- .find_bubbles_in_tx_pair(left_SSids, right_SSids)
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
        spath <- spath(x, gene_id=gene_id)
        bubbles(spath)
    }
)

setMethod("bubbles", "IntegerList",
    function(x, gene_id=NA)
    {
        if (!identical(gene_id, NA))
            stop("the 'gene_id' arg is not supported ",
                 "when 'x' is an IntegerList")
        .extract_bubbles_from_spath(x)
    }
)

### TODO: Add "bubbles" methods for data.frame and DataFrame objects.

