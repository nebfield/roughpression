#!/usr/bin/env Rscript

set.seed(0451)

library("dada2")
library("phyloseq")

# continue pipeline https://benjjneb.github.io/dada2/tutorial_1_4.html

# parameters
# -----------------------------------------------------------------------------
# 1: Rdata file from dada2_oral
# 2: tab separated file of sample data 

args <- commandArgs(trailingOnly=TRUE)
load(args[[1]])

# assign taxonomy
download.file(
  "https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz",
  destfile = "silva_128_train.fa.gz"
)

taxa <- assignTaxonomy(seqtab_nochim, "silva_128_train.fa.gz", multithread=8)

# assign species
download.file(
  "https://zenodo.org/record/824551/files/silva_species_assignment_v128.fa.gz",
  destfile = "silva_species_assignment_v128.fa.gz"
)

taxa_plus <-
  addSpecies(taxa,
             "silva_species_assignment_v128.fa.gz",
             verbose = FALSE)

# add labels
samdf <-
  data.frame(read.table(
    file = args[[2]],
    sep = "\t",
    header = TRUE,
    na.strings = ""
  ),
  stringsAsFactors = FALSE)
samdf$id <- paste0("sample", samdf$id) 
rownames(samdf) <- samdf$id
samdf$cohort <- factor(samdf$cohort, labels = c("control", "depressed"))
samdf$sex <- factor(samdf$sex, labels = c("Male", "Female"))
samdf$SMOKING <- factor(samdf$SMOKING, labels = c("Past", "Daily", "Occasionally", "Never"))
samdf$id <- NULL

oral <- phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa_plus)
)

oral <- phyloseq::subset_taxa(oral, !is.na(Phylum))

saveRDS(object = oral, file = "oral.rds")