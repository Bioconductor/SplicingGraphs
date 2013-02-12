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

library(SplicingGraphs)
TSPC_utils_path <- system.file("scripts", "TSPC-utils.R",
                               package="SplicingGraphs", mustWork=TRUE)
source(TSPC_utils_path)

### Make a TSPC splicing graph data frame and save it in the current working
### directory.
makeAndSaveTSPCsgdf <- function(subdir_path)
{
    subdir_basename <- basename(subdir_path)
    objname <- paste0(subdir_basename, "sgdf")
    filename <- paste0(objname, ".rda")
    sgdf <- makeTSPCsgdf(subdir_path)
    message("Saving ", objname, " to ", filename, " ... ", appendLF=FALSE)
    assign(objname, sgdf, envir=.GlobalEnv)
    save(list=objname, file=filename, envir=.GlobalEnv)
    message("OK")
}

makeAndSaveAllTSPCsgdfs <- function(subdir_paths)
{
    for (subdir_path in subdir_paths)
        makeAndSaveTSPCsgdf(subdir_path)
}

### Run this to make and save all the TSPC splicing graph data frames:
#makeAndSaveAllTSPCsgdfs(subdir_paths)

