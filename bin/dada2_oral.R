#!/usr/bin/env Rscript
set.seed(0451)

library("dada2"); packageVersion("dada2")
library("ShortRead")
library("ggplot2")

fnFs <- sort(list.files(getwd(), pattern="_R1.fastq.gz"))
fnRs <- sort(list.files(getwd(), pattern="_R2.fastq.gz"))
samp_names <- paste0("sample", sapply(strsplit(fnFs, "_"), `[`, 1))

filtFs <- paste0(samp_names, "_F_filt.fastq.gz")
filtRs <- paste0(samp_names, "_R_filt.fastq.gz")

out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240, 225),
                            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE)
#
errF <- dada2::learnErrors(filtFs, multithread=TRUE)
errR <- dada2::learnErrors(filtRs, multithread=TRUE)
plot_errF <- dada2::plotErrors(errF, nominalQ = TRUE)
ggsave("forward-error.pdf", plot_errF, device = "pdf")

derepFs <- dada2::derepFastq(filtFs, verbose=TRUE)
derepRs <- dada2::derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- samp_names
names(derepRs) <- samp_names

dadaFs <- dada2::dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(derepRs, err=errR, multithread=TRUE)

mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- dada2::makeSequenceTable(mergers)

# homings says ~441 basepairs
# distribution shows a couple thousand reads below 420bp, "cut a band"
distribution <- tibble::as_tibble(table(nchar(colnames(seqtab))))
colnames(distribution) <- c("length", "frequency")
distribution$reference <- 441
write.table(distribution, file = "seqlen-distribution.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
seqtab_trim <- seqtab[,nchar(colnames(seqtab)) %in% seq(423,431)]

seqtab_nochim <-
  dada2::removeBimeraDenovo(seqtab_trim,
                            method = "consensus",
                            multithread = TRUE,
                            verbose = TRUE)
save(seqtab_nochim, file = "seqtab.rda")

# sanity check
getN <- function(x) sum(dada2::getUniques(x))
track <-
  cbind(
    out,
    sapply(dadaFs, getN),
    sapply(mergers, getN),
    rowSums(seqtab),
    rowSums(seqtab_trim),
    rowSums(seqtab_nochim)
  )
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled",
                     "trimmed", "nonchim")
rownames(track) <- samp_names
write.table(track, file = "sanity_check.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
save.image(file = "dada2.RData")
