### =========================================================================
### SplicingGraphs objects
### -------------------------------------------------------------------------


### We deliberately choose to not extend GRangesList to make SplicingGraphs
### objects read-only and with a very restricted API (opaque objects).
setClass("SplicingGraphs",
    representation(
        tx="GRangesList"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Basic accessors.
###
### We support only a very small subset of getters from the GRangesList API.
###

setMethod("length", "SplicingGraphs", function(x) length(x@tx))

setMethod("names", "SplicingGraphs", function(x) names(x@tx))

setMethod("elementLengths", "SplicingGraphs", function(x) elementLengths(x@tx))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

setMethod("show", "SplicingGraphs",
    function(object)
    {
        ntx <- length(object)
        object_names <- names(object)
        if (is.null(object_names)) {
            ngenes <- ifelse(ntx == 0L, 0L, 1L)
        } else {
            ngenes <- length(unique(object_names))
        }
        cat(class(object), " object with ", ngenes, " gene(s) ",
            "and ", ntx, " transcript(s)\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SplicingGraphs() constructor
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

    ## Set spath metadata col.
    if ("spath" %in% colnames(exbytx_mcols))
        stop("'exbytx' already has metadata column spath")
    if (on.minus.strand) {
        spath <- rbind(SSids$end_SSid, SSids$start_SSid)
    } else {
        spath <- rbind(SSids$start_SSid, SSids$end_SSid)
    }
    spath_partitioning <- PartitioningByEnd(end(PartitioningByEnd(exbytx)) * 2L)
    names(spath_partitioning) <- tx_id
    spath <- splitAsList(as.vector(spath), spath_partitioning)
    exbytx_mcols$spath <- spath

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
###       are considered to belong to the same gene and the SplicingGraphs
###       object returned by SplicingGraphs() will be unnamed.
###   (b) A list of integer or character vectors, or an IntegerList, or a
###       CharacterList object, of length the number of genes to process,
###       and where 'grouping[[i]]' is a vector of valid subscripts in 'exbytx'
###       pointing to all the transcripts of the i-th gene.
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
.make_SplicingGraphs_from_GRangesList <- function(exbytx, grouping=NULL,
                                                  check.introns=TRUE)
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

### TODO: Make this a generic function (and rename 'exbytx' arg -> 'x').
SplicingGraphs <- function(exbytx, grouping=NULL, check.introns=TRUE)
{
    ans_tx <- .make_SplicingGraphs_from_GRangesList(exbytx, grouping=grouping,
                                                    check.introns=check.introns)
    new("SplicingGraphs", tx=ans_tx)
}

