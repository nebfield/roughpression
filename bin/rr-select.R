#!/usr/bin/env Rscript

set.seed(0451)

# 1: phyloseq object
# 2: feature selection text file 
args <- commandArgs(trailingOnly=TRUE)

ps <- readRDS(args[[1]])

# java uses zero-based indexing for arrays
# R uses one-based indexing for arrays
idx <-
  as.numeric(read.table(
    args[[2]],
    sep = ",",
    header = FALSE
  ))
idx_fixed <- idx + 1
 
# trim <- phyloseq::otu_table(ps)[, idx_fixed]
# trim_ps <- phyloseq::prune_taxa(colnames(trim), ps)
# saveRDS(object = trim_ps, file = "trimmed.rds")
