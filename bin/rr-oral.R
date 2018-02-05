#!/usr/bin/env Rscript

set.seed(0451)
args = commandArgs(trailingOnly=TRUE)

oral <- readRDS(args[[1]])
oral <- phyloseq::transform_sample_counts(oral, function(x) x / sum(x))
oral <- phyloseq::subset_samples(oral, SMOKING == "Never")
oral_y <- as.numeric(phyloseq::sample_data(oral)$cohort)

# objects: rows, attributes: columns
ra_annot <- cbind(data.frame(phyloseq::otu_table(oral), Class = oral_y, check.names = FALSE))
write.table(ra_annot, file = "oral.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
