AScode2desc <- NULL

.onLoad <- function(libname, pkgname)
{
    filepath <- system.file("extdata", "ASpatterns.txt",
                            package=pkgname, lib.loc=libname,
                            mustWork=TRUE)
    ASpatterns <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE)
    tmp <- ASpatterns[ , "Description"]
    names(tmp) <- ASpatterns[ , "AScode"]
    AScode2desc <<- tmp
}

