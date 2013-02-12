### The TSPC data is on the rhinos.

TSPC_path <- "/shared/labs/mmcintos/proj/Solexa/TSPC"
TSPC_subdirs <- c("BAI1", "CYB561", "DAPL1", "ITGB8", "KIAA0319L",
                  "LGSN", "MKRN3", "MUC16", "ST14", "TREM2")

#stopifnot(identical(file.path(TSPC_path, TSPC_subdirs), setdiff(list.dirs(TSPC_path), TSPC_path)))

### No BAM files for KIAA0319L and TREM2.
### Transcripts T-4 and T-5 in MUC16 both have their 2nd exon included in
### their 3rd exon ==> splicing graph theory doesn't apply.
exclude_subdirs <- c("KIAA0319L", "TREM2", "MUC16")

keep_subdirs <- setdiff(TSPC_subdirs, exclude_subdirs)
subdir_paths <- file.path(TSPC_path, keep_subdirs)

### Sanity check.

library(Rsamtools)

flag0 <- scanBamFlag(#isProperPair=TRUE,
                     isNotPrimaryRead=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isDuplicate=FALSE)
param0 <- ScanBamParam(flag=flag0, what="flag")

library(SplicingGraphs)
TSPC_utils_path <- system.file("scripts", "TSPC-utils.R",
                               package="SplicingGraphs", mustWork=TRUE)
source(TSPC_utils_path)

makeSgdf <- function(subdir_path)
{
    subdir_basename <- basename(subdir_path)
    filenames <- list.files(subdir_path)
    filenames_nchar <- nchar(filenames)

    ## Load the gene model.
    suffixes <- substr(filenames, filenames_nchar-10L, filenames_nchar)
    models_filename <- filenames[suffixes == "_models.txt"]
    models_path <- file.path(subdir_path, models_filename)
    cat("Reading ", models_path, " ...", sep="")
    ex_by_tx <- loadModels(models_path)
    cat(" OK\n")

    ## Compute the splicing graph.
    ex_by_tx2 <- splicingGraphs(ex_by_tx)

    ## Load each BAM file
    suffixes <- substr(filenames, filenames_nchar-3L, filenames_nchar)
    bam_filenames <- filenames[suffixes == ".bam"]
    prefixes <- substr(bam_filenames, 1L, nchar(subdir_basename)+1L)
    stopifnot(all(prefixes == paste0(subdir_basename, "-")))
    sample_names <- substr(bam_filenames, nchar(prefixes)+1L,
                                          nchar(bam_filenames)-4L)
    nbam <- length(bam_filenames)
    cat(nbam, " BAM files to process\n", sep="")
    nread <- integer(nbam)
    X <- seq_len(nbam)
    names(X) <- sample_names
    nhits <- sapply(X, function(i) {
        bam_filename <- bam_filenames[i]
        sample_name <- sample_names[i]
        bam_filepath <- file.path(subdir_path, bam_filename)
        gal <- readGappedAlignments(bam_filepath, use.names=TRUE,
                                    param=param0)
        is_paired <- bamFlagTest(mcols(gal)$flag, "isPaired")
        if (!any(is_paired)) {
            nread[i] <- length(gal)
            grl <- grglist(gal, order.as.in.query=TRUE)
        } else {
            stopifnot(all(is_paired))
            galp <- readGappedAlignmentPairs(bam_filepath, use.names=TRUE,
                                             param=param0)
            nread[i] <- length(galp)
            grl <- grglist(galp, order.as.in.query=TRUE)
        }
        sgdf <- makeSgdfWithHits(grl, ex_by_tx2)
        sgdf[ , "nhits"]
    })
    DataFrame(Sgdf(ex_by_tx2), nhits)
}

for (subdir_path in subdir_paths) {
    makeSgdf(subdir_path)
}

