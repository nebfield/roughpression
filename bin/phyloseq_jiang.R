#!/usr/bin/env Rscript

library("dada2")
library("phyloseq")

# parameters
# -----------------------------------------------------------------------------
# 1: Rdata file from dada2_jiang

args <- commandArgs(trailingOnly=TRUE)
load(args[[1]])

# continue pipeline https://benjjneb.github.io/dada2/tutorial_1_4.html

seqtab <- makeSequenceTable(dadaFs)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(478,489)] # cut a band

seqtab.nochim <-
  removeBimeraDenovo(seqtab2,
                     method = "consensus",
                     multithread = TRUE,
                     verbose = TRUE)

# sanity check
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- rownames(seqtab2)
head(track)
write.table(track, file = "sanity-check.txt", sep = "\t", quote = FALSE)

# assign taxonomy
download.file(
  "https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz",
  destfile = "silva_128_train.fa.gz"
)

taxa <- assignTaxonomy(seqtab.nochim, "silva_128_train.fa.gz", multithread=TRUE)

samdf <-
  data.frame(cohort = sapply(substr(rownames(seqtab.nochim), 1, 1), function(x)
    if (x == "C") {
      return("Control")
    } else if (x == "M") {
      return("MDD")
    } else if (x == "R") {
      return("Remission")
    }))
rownames(samdf) <- rownames(seqtab.nochim)

gut <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)

saveRDS(object = gut, file = "gut.rds")
