#!/usr/bin/env Rscript

library(dada2); packageVersion("dada2")

# args: list of filenames
args <- commandArgs(trailingOnly=TRUE)

fnFs <- args
# Extract sample names, assuming filenames have format: XXX_SAMPLENAME.fastq.gz
sampleNames <- gsub("\\.fastq\\.gz", "", gsub("^.*_", "", fnFs))
filtFs <- file.path(paste0(sampleNames, "_filt.fastq.gz"))

# V1 - V3 region per paper
# should be ~600 nucleotides long
out <-
  dada2::filterAndTrim(
    fwd = fnFs,
    filt = filtFs,
    maxN = 0,
    maxEE = 2,
    maxLen = 650,
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  ) 

errF <- dada2::learnErrors(filtFs, multithread=TRUE)
derepFs <- dada2::derepFastq(filtFs, verbose=FALSE)
names(derepFs) <- sampleNames
# extra args for 454 data
dadaFs <- dada2::dada(derepFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
save.image(file = "dada2.Rdata")