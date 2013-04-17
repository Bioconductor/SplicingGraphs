### =========================================================================
### makeTSPCsg.R -- Script for making the SplicingGraphs object for the TSPC
###                 project
### -------------------------------------------------------------------------

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

### Make the SplicingGraphs object.
TSPCsg <- make_TSPC_SplicinGraphs(subdir_paths)

### Make the BAM status matrix.
### BAM status:
###   ".": BAM file doesn't exist;
###   "0": file is empty (no alignments);
###   "s": single-end;
###   "p": paired-end;
###   "m": mixed single-/paired-end.
sample_names <- get_TSPC_sample_names(subdir_paths)
bam_status_matrix <- make_TSPC_bam_status_matrix(subdir_paths, sample_names)
dim(bam_status_matrix)
bam_status_matrix[ , 1:8]
### A close look at the matrix reveals that:
###   - 9 TSPC genes, only 7 with BAM files: BAI1, CYB561, DAPL1, ITGB8, LGSN,
###     MKRN3, and ST14. No BAM files for KIAA0319L and TREM2.
###   - 54 TSPC samples: 42 have single-end reads, 12 have paired-end reads.

### Make the "BAM gap rate" matrix.
bam_gaprate_matrix <- make_TSPC_bam_gaprate_matrix(subdir_paths, sample_names)m_gaprate_matrix[ , 1:8]

### Assign the reads to the SplicingGraphs object (this takes about 4 min on
### rhino02).
TSPCsg <- assign_TSPC_reads(TSPCsg, subdir_paths)

### Serialized the SplicingGraphs object.
save(TSPCsg, file="TSPCsg.rda", compress="xz")

