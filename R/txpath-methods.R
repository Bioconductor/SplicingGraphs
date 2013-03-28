### =========================================================================
### "txpath" (and related) methods
### -------------------------------------------------------------------------


get_sgnodes_from_txpath <- function(txpath)
{
    SSids <- unique(unlist(txpath, use.names=FALSE))
    c("R", sort(SSids), "L")
}

make_matrix_from_txpath <- function(txpath)
{
    sgnodes <- get_sgnodes_from_txpath(txpath)
    ans_nrow <- length(txpath)
    ans_ncol <- length(sgnodes)
    ans_dimnames <- list(names(txpath), sgnodes)
    ans <- matrix(FALSE , nrow=ans_nrow, ncol=ans_ncol, dimnames=ans_dimnames)
    ans[ , 1L] <- ans[ , ans_ncol] <- TRUE
    i <- cbind(rep.int(seq_along(txpath), elementLengths(txpath)),
               unlist(txpath, use.names=FALSE) + 1L)
    ans[i] <- TRUE
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### txpath() accessor
###
### Gets the splicing paths.
### Returns them in a named IntegerList object (actually a
### CompressedIntegerList instance) with 1 list element per transcript.
### Each list element 'txpath[[i]]' contains the Splicing Site ids for the
### i-th transcript.
###

setGeneric("txpath", signature="x",
    function(x, as.matrix=FALSE) standardGeneric("txpath")
)

setMethod("txpath", "GRangesList",
    function(x, as.matrix=FALSE)
    {
        if (!identical(as.matrix, FALSE))
            stop("the 'as.matrix' arg is not supported ",
                 "when 'x' is a GRangesList object")
        ans <- mcols(x)[ , "txpath"]
        x_names <- names(x)
        if (!is.null(x_names))
            names(ans) <- x_names
        ans
    }
)

setMethod("txpath", "SplicingGraphs",
    function(x, as.matrix=FALSE)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        if (!isTRUEorFALSE(as.matrix))
            stop("'as.matrix' must be TRUE or FALSE")
        ## Same as 'x[[1]]' but should be faster.
        x_1 <- unlist(x, use.names=FALSE)
        ans <- txpath(x_1)
        if (!as.matrix)
            return(ans)
        make_matrix_from_txpath(ans)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UATXHcount() accessor
###

setGeneric("UATXHcount", signature="x",
    function(x) standardGeneric("UATXHcount")
)

### Should return an integer vector or a NULL.
setMethod("UATXHcount", "SplicingGraphs",
    function(x)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        mcols(unlist(x, use.names=FALSE))[["UATXHcount"]]
    }
)

