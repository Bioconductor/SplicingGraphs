### =========================================================================
### SplicingGraphs objects
### -------------------------------------------------------------------------


### A simpler design would be to define only 1 class, the SplicingGraphs
### class, and to define it the way the .SplicingGraphGenes class below is
### defined. The problem with this is that the SplicingGraphs class then
### inherits a very rich API (Vector + List) but many operations (like c()
### or relist()) are broken, unless we implement specific methods for them.
### But: (a) that's a lot of work (the API is huge), and (b) we don't need
### those operations in the first place. All we need are: length(), names(),
### [, [[, elementLengths(), and unlist(). Hence the 2 class definitions
### below. The .SplicingGraphGenes class is an internal class that is not
### intended to be exposed to the user.

setClass(".SplicingGraphGenes",
    contains="CompressedList",
    representation(
        unlistData="GRangesList",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GRangesList"
    )
)

setClass("SplicingGraphs",
    representation(
        genes=".SplicingGraphGenes"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###
### A valid SplicingGraphs object is an object with a valid genes slot. So
### it's enough to implement a validity method for .SplicingGraphGenes
### objects and to validate a SplicingGraphs object 'x' with:
###
###   validObject(x, complete=TRUE)
###

.valid.SplicingGraphGenes.names <- function(x)
{
    x_names <- names(x)
    if (is.null(x_names)) {
        if (length(x) == 1L)
            return(NULL)
        return("'x' must have names")
    }
    if (anyDuplicated(x_names))
        return("'x' has duplicated names")
    NULL
}

.valid.SplicingGraphGenes.unlistData <- function(x)
{
    x_unlistData <- x@unlistData
    if (!is.null(names(x_unlistData)))
        return("'x@unlistData' must be unnamed")
    NULL
}

.valid.SplicingGraphGenes <- function(x)
{
    c(.valid.SplicingGraphGenes.names(x),
      .valid.SplicingGraphGenes.unlistData(x))
}

setValidity2(".SplicingGraphGenes", .valid.SplicingGraphGenes)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Restricted SplicingGraphs API.
###

setMethod("length", "SplicingGraphs", function(x) length(x@genes))

setMethod("names", "SplicingGraphs", function(x) names(x@genes))

setMethod("[", "SplicingGraphs",
    function(x, i, j, ... , drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        x@genes <- x@genes[i, drop=drop]
        x
    }
)

setMethod("[[", "SplicingGraphs",
    function (x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        x@genes[[i]]

    }
)

setMethod("elementLengths", "SplicingGraphs",
    function(x) elementLengths(x@genes)
)

setMethod("unlist", "SplicingGraphs",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        unlist(x@genes, recursive=recursive, use.names=use.names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

setMethod("show", "SplicingGraphs",
    function(object)
    {
        ngene <- length(object)
        ntx <- length(unlist(object, use.names=FALSE))
        cat(class(object), " object with ", ngene, " gene(s) ",
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

### 'gene' must be a GRangesList object containing the exons of a *single*
### gene grouped by transcripts. More precisely, each top-level element
### in 'gene' contains the genomic ranges of the exons for a particular
### transcript of the gene.
.setSplicingGraphInfo <- function(gene, check.introns=TRUE)
{
    if (!is(gene, "GRangesList"))
        stop("'gene' must be a GRangesList object")
    if (!isTRUEorFALSE(check.introns))
        stop("'check.introns' must be TRUE or FALSE")
    exons <- gene@unlistData
    exons_strand <- strand(exons)
    if (nrun(seqnames(exons)) != 1L || nrun(exons_strand) != 1L)
        stop("all the exons in the gene must be on the same ",
             "reference sequence and strand")
    on.minus.strand <- runValue(exons_strand)[1L] == "-"
    if (check.introns) {
        ## We check that, within each transcript, exons are ordered from 5'
        ## to 3' with gaps of at least 1 nucleotide between them.
        ranges_by_tx <- ranges(gene)
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
        stop("'unlist(gene)' already has metadata columns: ", in_1_string)
    }
    mcols(exons) <- cbind(mcols(exons), DataFrame(SSids))
    gene@unlistData <- exons
    gene_mcols <- mcols(gene)

    ## Set tx_id metadata col.
    if ("tx_id" %in% colnames(gene_mcols))
        stop("'gene' already has metadata column tx_id")
    tx_id <- names(gene)
    if (!is.null(tx_id))
        gene_mcols$tx_id <- tx_id

    ## Set txpaths metadata col.
    if ("txpaths" %in% colnames(gene_mcols))
        stop("'gene' already has metadata column txpaths")
    if (on.minus.strand) {
        txpaths <- rbind(SSids$end_SSid, SSids$start_SSid)
    } else {
        txpaths <- rbind(SSids$start_SSid, SSids$end_SSid)
    }
    txpaths_partitioning_end <- end(PartitioningByEnd(gene)) * 2L
    txpaths_partitioning <- PartitioningByEnd(txpaths_partitioning_end)
    names(txpaths_partitioning) <- tx_id
    txpaths <- splitAsList(as.vector(txpaths), txpaths_partitioning)
    gene_mcols$txpaths <- txpaths

    mcols(gene) <- gene_mcols
    gene
}

### 'x' must be a GRangesList object containing the exons of one or more
### genes grouped by transcripts. More precisely, each top-level element
### in 'x' contains the genomic ranges of the exons for a particular
### transcript. Typically 'x' will be obtained from a TranscriptDb object
### 'txdb' with 'exonsBy(txdb, by="tx")'.
### 'grouping' is an optional object that represents the grouping by gene of
### the top-level elements (i.e. transcripts) in 'x'. It can be either:
###   (a) Missing (i.e. NULL). In that case, all the transcripts in 'x'
###       are considered to belong to the same gene and the SplicingGraphs
###       object returned by SplicingGraphs() will be unnamed.
###   (b) A list of integer or character vectors, or an IntegerList, or a
###       CharacterList object, of length the number of genes to process,
###       and where 'grouping[[i]]' is a vector of valid subscripts in 'x'
###       pointing to all the transcripts of the i-th gene.
###   (c) A factor, character vector, or integer vector, of length 'x'
###       with 1 level per gene.
###   (d) A named GRangesList object containing transcripts grouped by genes
###       i.e. each top-level element in 'grouping' contains the genomic ranges
###       of the transcripts for a particular gene. In that case, the grouping
###       is inferred from the tx_id (or alternatively tx_name) metadata
###       column of 'unlist(grouping)' and all the values in that column must
###       be in 'names(x)'.
###       If 'x' was obtained with 'exonsBy(txdb, by="tx")', then the
###       GRangesList object used for grouping would typically be obtained with
###       'transcriptsBy(txdb, by="gene")'.
###   (e) A data.frame or DataFrame with 2 character vector columns: a
###       gene_id column (factor, character vector, or integer vector),
###       and a tx_id (or alternatively tx_name) column. In that case, 'x'
###       must be named and all the values in the tx_id (or tx_name) column
###       must be in 'names(x)'.

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
.normargGrouping <- function(grouping, x)
{
    ## (b)
    if (is.list(grouping) || is(grouping, "IntegerList")
     || is(grouping, "CharacterList")) {
        return(.checkOrMakeUniqueGroupingNames(grouping))
    }
    ## (c)
    if (is.factor(grouping) || is.character(grouping) || is.integer(grouping)) {
        if (length(grouping) != length(x))
            stop("when 'grouping' is a factor, character vector, or integer ",
                 "vector, it must have the same length as 'x'")
        return(split(seq_along(x), grouping))
    }
    x_names <- names(x)
    ## (d)
    if (is(grouping, "GRangesList")) {
        if (is.null(x_names))
            stop("when 'grouping' is a GRangesList, 'x' must be named ",
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
            m <- match(mcols[[idx]], x_names)
            if (any(is.na(m)))
                next
            return(splitAsList(m, PartitioningByEnd(grouping)))
        }
        stop("'unlist(grouping)' has no tx_id or tx_name column, ",
             "or they contain values that are not in 'names(x)'")
    }
    ## (e)
    if (is.data.frame(grouping) || is(grouping, "DataFrame")) {
        if (is.null(x_names))
            stop("when 'grouping' is a data.frame or a DataFrame, 'x' ",
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
            m <- match(grouping[[idx]], x_names)
            if (any(is.na(m)))
                next
            return(split(m, gene_id))
        }
        stop("'grouping' has no tx_id or tx_name column, ",
             "or they contain values that are not in 'names(x)'")
    }
    stop("invalid 'grouping'")

}

### TODO: Improve handling of invalid genes i.e. provide more details about
### which genes were considered invalid and why.
.make_unlisted_genes_from_GRangesList <- function(x, grouping=NULL,
                                                     min.ntx=2, max.ntx=NA,
                                                     check.introns=TRUE)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    if (is.null(grouping)) {
        if (!identical(min.ntx, 2) || !identical(max.ntx, NA))
            stop("the 'min.ntx' and 'max.ntx' args are not supported ",
                 "when 'grouping' is not supplied or NULL")
        ans <- .setSplicingGraphInfo(x, check.introns=check.introns)
        names(ans) <- NULL
        return(ans)
    }

    ## Check 'min.ntx'.
    if (!isSingleNumber(min.ntx))
        stop("'min.ntx' must be a single number")
    if (!is.integer(min.ntx))
        min.ntx <- as.integer(min.ntx)
    if (min.ntx < 1L)
        stop("'min.ntx' must be >= 1")

    ## Check 'max.ntx'.
    if (!isSingleNumberOrNA(max.ntx))
        stop("'max.ntx' must be a single number or NA")
    if (!is.integer(max.ntx))
        max.ntx <- as.integer(max.ntx)
    if (!is.na(max.ntx) && max.ntx < min.ntx)
        stop("'max.ntx' must be >= 'min.ntx'")

    grouping <- .normargGrouping(grouping, x)

    ## Keep genes with nb of transcripts >= min.ntx and <= max.ntx.
    grouping_eltlen <- elementLengths(grouping)
    keep <- grouping_eltlen >= min.ntx
    if (!is.na(max.ntx))
        keep <- keep & grouping_eltlen <= max.ntx
    grouping <- grouping[keep]

    ## Main loop.
    ans <- lapply(seq_along(grouping),
                  function(i) {
                      ii <- grouping[[i]]
                      gene <- x[ii]
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

.new_SplicingGraphGenes <- function(unlisted_genes)
{
    unlisted_names <- names(unlisted_genes)
    if (is.null(unlisted_names)) {
        ans_partitioning <- PartitioningByEnd(length(unlisted_genes))
    } else {
        names(unlisted_genes) <- NULL
        ans_end <- end(Rle(unlisted_names))
        ans_names <- unlisted_names[ans_end]
        ans_partitioning <- PartitioningByEnd(ans_end, names=ans_names)
    }
    IRanges:::newCompressedList0(".SplicingGraphGenes",
                                 unlisted_genes, ans_partitioning)
}

setGeneric("SplicingGraphs", signature="x",
    function(x, grouping=NULL, min.ntx=2, max.ntx=NA, check.introns=TRUE)
        standardGeneric("SplicingGraphs")
)

setMethod("SplicingGraphs", "GRangesList",
    function(x, grouping=NULL, min.ntx=2, max.ntx=NA, check.introns=TRUE)
    {
        unlisted_genes <- .make_unlisted_genes_from_GRangesList(x,
                            grouping=grouping, min.ntx=min.ntx, max.ntx=max.ntx,
                            check.introns=check.introns)
        ans_genes <- .new_SplicingGraphGenes(unlisted_genes)
        new("SplicingGraphs", genes=ans_genes)
    }
)

setMethod("SplicingGraphs", "TranscriptDb",
    function(x, grouping=NULL, min.ntx=2, max.ntx=NA, check.introns=TRUE)
    {
        if (!is.null(grouping))
            stop("the 'grouping' arg is not supported ",
                 "when 'x' is a TranscriptDb object")
        ex_by_tx <- exonsBy(x, by="tx", use.names=TRUE)
        tx_by_gn <- transcriptsBy(x, by="gene")
        SplicingGraphs(ex_by_tx, tx_by_gn, min.ntx=min.ntx, max.ntx=max.ntx,
                       check.introns=check.introns)
    }
)

