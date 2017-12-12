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

# seqtab <- dadaFs
# # Inspect distribution of sequence lengths
# table(nchar(dada2::getSequences(seqtab)))
# seqtab_nochim <-
#   dada2::removeBimeraDenovo(seqtab,
#                             method = "consensus",
#                             multithread = TRUE,
#                             verbose = TRUE)

# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab_nochim))
# colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
# rownames(track) <- sampleNames
# head(track)
# 
# download.file(
#   url = "https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz",
#   destfile = "silva_128_train.fa.gz"
# )
# taxa <- assignTaxonomy(seqtab.nochim, "silva_128_train.fa.gz", multithread=TRUE)
# 
# saveRDS(object = seqtab_nochim, file = "jiang_seqtab.rds")
# saveRDS(object = taxa, file = "jiang_taxa.rds")