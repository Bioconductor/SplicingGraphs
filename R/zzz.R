ASCODE2DESC <- NULL

.onLoad <- function(libname, pkgname)
{
    ## We need to fix the prototypes of the GeneModel and SplicingGraphs
    ## classes. Without these fixes, 'new("GeneModel")' and
    ## 'new("SplicingGraphs")' would return invalid objects e.g.
    ##
    ##     validObject(new("SplicingGraphs"), complete=TRUE)
    ##
    ## would fail.
    sg0 <- emptySplicingGraphs()
    IRanges:::setPrototypeFromObject("GeneModel",
                                     sg0@genes,
                                     where=asNamespace(pkgname))
    IRanges:::setPrototypeFromObject("SplicingGraphs",
                                     sg0,
                                     where=asNamespace(pkgname))

    ## Set the ASCODE2DESC global constant based on the content of the
    ## extdata/ASpatterns.txt file.
    filepath <- system.file("extdata", "ASpatterns.txt",
                            package=pkgname, lib.loc=libname,
                            mustWork=TRUE)
    ASpatterns <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE)
    tmp <- ASpatterns[ , "Description"]
    names(tmp) <- ASpatterns[ , "AScode"]
    ASCODE2DESC <<- tmp
}

