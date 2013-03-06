
.make_spath_matrix_from_spath <- function(spath)
{
    nodes <- c("R", sort(unique(unlist(spath))), "L")
}

findBubbles <- function(x, gene_id=NA)
{
    spath <- spath(sg, gene_id=gene_id)
    spath_mat <- .make_spath_matrix_from_spath(spath)
    spath_mat
}

