### The TSPC data is on the rhinos.

TSPC_path <- "/shared/labs/mmcintos/proj/Solexa/TSPC"
TSPC_subdirs <- c("BAI1", "CYB561", "DAPL1", "ITGB8", "KIAA0319L",
                  "LGSN", "MKRN3", "MUC16", "ST14", "TREM2")

#stopifnot(identical(file.path(TSPC_path, TSPC_subdirs), setdiff(list.dirs(TSPC_path), TSPC_path)))

### Transcripts T-4 and T-5 in MUC16 both have their 2nd exon included in
### their 3rd exon ==> splicing graph theory doesn't apply.
exclude_subdirs <- "MUC16"

keep_subdirs <- setdiff(TSPC_subdirs, exclude_subdirs)
subdir_paths <- file.path(TSPC_path, keep_subdirs)

library(SplicingGraphs)
TSPC_utils_path <- system.file("scripts", "TSPC-utils.R",
                               package="SplicingGraphs", mustWork=TRUE)
source(TSPC_utils_path)

TSPCsg <- make_TSPC_SplicinGraphs(subdir_paths)
TSPCsg <- assign_TSPC_reads(TSPCsg, subdir_paths)
save(TSPCsg, file="TSPCsg.rda", compress="xz")

