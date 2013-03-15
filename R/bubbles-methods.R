### =========================================================================
### "bubbles" methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_sgnodetypes_from_txpathmat()
###

### Supported types are 1 (exon start), 2 (exon end), 0 (R or L nodes).
### Returns an integer vector of length the nb of cols in 'txpathmat' and
### named with 'txpathmat' colnames.
### TODO: (a) Add an sgnodetypes() accessor to sgedges-methods.R, (b) move
### the .get_sgnodetypes_from_txpathmat() helper to that file, and (c) use
### it internally in sgnodetypes().
.get_sgnodetypes_from_txpathmat <- function(txpathmat, check.matrix=FALSE)
{
    ans <- integer(ncol(txpathmat))
    names(ans) <- colnames(txpathmat)
    for (i in seq_len(nrow(txpathmat))) {
        idx <- which(txpathmat[i, , drop=FALSE])
        if (length(idx) <= 2L)
            next()
        idx <- idx[-c(1L, length(idx))]
        exon_start_idx <- idx[c(TRUE, FALSE)]
        exon_end_idx <- idx[c(FALSE, TRUE)]
        if (check.matrix) {
            if (any(ans[exon_start_idx] == 2L) || any(ans[exon_end_idx] == 1L))
                stop("invalid matrix of transcript paths: ",
                     "some columns in 'txpathmat' seem to correspond ",
                     "at the same time to an exon start and an exon end")
        }
        ans[exon_start_idx] <- 1L
        ans[exon_end_idx] <- 2L
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_bubbles_from_txpathmat()
###

### We must have 1 <= i, j <= ncol(txpathmat), and j - i >= 2. This is not
### checked!
.is_bubble <- function(txpathmat, i, j)
{
    txbase <- txpathmat[, i] & txpathmat[, j]
    if (sum(txbase) <= 1L)
        return(FALSE)
    for (k in (i+1L):(j-1L))
        if (all(txpathmat[ , k] >= txbase))
            return(FALSE)
    TRUE
}

### A fancy alternative to .is_bubble() that tries to avoid the for loop.
### Surprisingly, it turned out to be slower than .is_bubble() in practise.
.is_bubble2 <- function(txpathmat, i, j)
{
    txbase <- txpathmat[, i] & txpathmat[, j]
    if (sum(txbase) <= 1L)
        return(FALSE)
    all(colSums(txpathmat[ , (i+1L):(j-1L), drop=FALSE] < txbase) >= 1)
}

### Assumes 'm' to be a logical matrix. Not checked.
.make_bit_strings_from_logical_matrix <- function(m)
{
    apply(m, 1L, function(x) paste(as.integer(x), collapse=""))
}

### Assumes 's' to be a vector of equal-length strings made of 0s and 1s. Not
### checked.
.make_logical_matrix_from_bit_strings <- function(s)
{
    all_letters <- unlist(strsplit(s, NULL, fixed=TRUE), use.names=FALSE)
    matrix(as.logical(as.integer(all_letters)), nrow=length(s), byrow=TRUE,
           dimnames=list(names(s), NULL))
}

### Returns a DataFrame with 1 row per variant and the following cols:
### partition <CharacterList>, path <CharacterList>, and code <character>.
.get_bubble_variants <- function(txpathmat, sgnodetypes, i, j)
{
    txbase <- txpathmat[, i] & txpathmat[, j]
    bubble_submat <- txpathmat[txbase, (i+1L):(j-1L), drop=FALSE]

    ## Remove cols with FALSEs only.
    bubble_submat <- bubble_submat[ , colSums(bubble_submat) != 0L, drop=FALSE]
    bubble_submat_rownames <- rownames(bubble_submat)
    bubble_submat_colnames <- colnames(bubble_submat)

    ## Compute variant paths (1 per row).
    ans_path <- CharacterList(
                    lapply(seq_len(nrow(bubble_submat)),
                           function(i)
                               bubble_submat_colnames[bubble_submat[i, ]])
                )

    ## Compute variant relative paths (1 per row).
    ans_rpath <- IntegerList(
                    lapply(seq_len(nrow(bubble_submat)),
                           function(i)
                               which(bubble_submat[i, ]))
                )

    ## Compute variant code (1 per row).
    ans_code <- sapply(seq_len(length(ans_path)),
                       function(k)
                       {
                           path <- ans_path[[k]]
                           path_len <- length(path)
                           if (path_len == 0L)
                               return("0")
                           types <- c("-", "^")[sgnodetypes[path]]
                           ## Sanity check would fail if 'sgnodetypes[path]'
                           ## contained 0s but this should never happen.
                           if (length(types) != path_len)
                               stop("some splicing sites in variant path ",
                                    "are of type 0")
                           if (i == 1L)
                               types[1L] <- "["
                           if (j == ncol(txpathmat))
                               types[length(types)] <- "]"
                           paste0(ans_rpath[[k]], types, collapse="")
                       })

    ## Order the variants by lexicographic order on their code.
    oo <- order(ans_code)
    bubble_submat_rownames <- bubble_submat_rownames[oo]
    ans_path <- ans_path[oo]
    #ans_rpath <- ans_rpath[oo]
    ans_code <- ans_code[oo]

    ## Identify unique variants.
    ans_code1 <- ans_code[-length(ans_code)]
    ans_code2 <- ans_code[-1L]
    is_not_dup <- c(TRUE, ans_code1 != ans_code2)
    ans_path <- ans_path[is_not_dup]
    #ans_rpath <- ans_rpath[is_not_dup]
    ans_code <- ans_code[is_not_dup]

    ## Compute variant partitions.
    ans_partition <- unname(splitAsList(bubble_submat_rownames,
                                        cumsum(is_not_dup)))

    ## Make and return the DataFrame.
    DataFrame(partition=ans_partition,
              path=ans_path,
              #rpath=ans_rpath,
              code=ans_code)
}

### Returns a DataFrame with 1 row per bubble and the following cols:
### source <character>, sink <character>, d <integer>, partitions
### <CharacterList>, paths <CharacterList>, AScode <character>, and
### description <character>.
#.extract_bubbles_from_txpathmat <- function(txpathmat, outdeg, indeg)
.extract_bubbles_from_txpathmat <- function(txpathmat)
{
    sgnodetypes <- .get_sgnodetypes_from_txpathmat(txpathmat)
    ans_source <- ans_sink <- ans_AScode <- character(0)
    ans_d <- integer(0)
    ans_partitions <- ans_paths <- CharacterList()
    ## Surprisingly, the outdeg/indeg trick turned out to not make any
    ## noticeable speed difference in practise.
    #ii0 <- which(outdeg >= 2L)
    #jj0 <- which(indeg >= 2L)
    #for (i in ii0) {
    #    jj <- jj0[jj0 >= i + 2L]
    #    for (j in jj) {
    ncol <- ncol(txpathmat)
    for (i in 1:(ncol-2L)) {
        for (j in (i+2L):ncol) {
            if (!.is_bubble(txpathmat, i, j))
                next
            bubble_variants <- .get_bubble_variants(txpathmat, sgnodetypes,
                                                    i, j)
            bubble_d <- nrow(bubble_variants)
            if (bubble_d <= 1L)
                next
            ans_source <- c(ans_source, colnames(txpathmat)[i])
            ans_sink <- c(ans_sink, colnames(txpathmat)[j])
            ans_d <- c(ans_d, bubble_d)
            ## Format the bubble partitions.
            bubble_partitions <- bubble_variants[ , "partition"]
            bubble_partitions <- sapply(bubble_partitions, paste, collapse=",")
            bubble_partitions <- paste0("{", bubble_partitions, "}")
            ans_partitions <- c(ans_partitions,
                                CharacterList(bubble_partitions))
            ## Format the bubble paths.
            bubble_paths <- bubble_variants[ , "path"]
            bubble_paths <- sapply(bubble_paths, paste, collapse=",")
            bubble_paths <- paste0("{", bubble_paths, "}")
            ans_paths <- c(ans_paths, CharacterList(bubble_paths))
            ## Format the bubble AScode.
            bubble_AScode <- paste(bubble_variants[ , "code"], collapse=",")
            ans_AScode <- c(ans_AScode, bubble_AScode)
        }
    }
    DataFrame(source=ans_source,
              sink=ans_sink,
              d=ans_d,
              partitions=ans_partitions,
              paths=ans_paths,
              AScode=ans_AScode,
              description=unname(ASCODE2DESC[ans_AScode]))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bubbles()
###

setGeneric("bubbles", signature="x",
    function(x) standardGeneric("bubbles")
)

setMethod("bubbles", "SplicingGraphs",
    function(x)
    {
        if (length(x) != 1L)
            stop("'x' must be a SplicingGraphs object of length 1")
        bubbles_cache <- x@.bubbles_cache
        ans <- try(get(names(x), envir=bubbles_cache, inherits=FALSE),
                   silent=TRUE)
        if (is(ans, "try-error")) {
            ans <- callNextMethod()
            assign(names(x), ans, envir=bubbles_cache)
        }
        ans
    }
)

setMethod("bubbles", "ANY",
    function(x)
    {
        txpaths <- txpaths(x)
        bubbles(txpaths)
    }
)

setMethod("bubbles", "IntegerList",
    function(x)
    {
        txpathmat <- make_matrix_from_txpaths(x)
        #outdeg <- outdeg(x)
        #indeg <- indeg(x)
        #.extract_bubbles_from_txpathmat(txpathmat, outdeg, indeg)
        .extract_bubbles_from_txpathmat(txpathmat)
    }
)

### TODO: Add "bubbles" methods for data.frame and DataFrame objects.

