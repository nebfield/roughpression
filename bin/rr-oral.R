#!/usr/bin/env Rscript

set.seed(0451)
args = commandArgs(trailingOnly=TRUE)

oral <- readRDS(args[[1]])
oral <- phyloseq::transform_sample_counts(oral, function(x) x / sum(x))
oral <- phyloseq::subset_samples(oral, SMOKING == "Never")
oral_y <- as.numeric(phyloseq::sample_data(oral)$cohort)

# csv file ---------------------------------------------------------------------
# objects: rows, attributes: columns
ra_annot <-
  cbind(data.frame(
    phyloseq::otu_table(oral),
    Class = oral_y,
    check.names = FALSE
  ))

write.table(
  ra_annot,
  file = "oral.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# arff file --------------------------------------------------------------------
short_names <-
  tibble::tibble(short = make.unique(stringr::str_sub(
    phyloseq::taxa_names(oral), start = 1, end = 5)), 
    long = phyloseq::taxa_names(oral))

ra_arff <-
  cbind(
    data.frame(
      phyloseq::otu_table(oral),
      Class = phyloseq::sample_data(oral)$cohort,
      check.names = FALSE
    )
  )

colnames(ra_arff) <- c(short_names$short, "Class")

write.table(
  short_names,
  file = "oral_short_names.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

foreign::write.arff(ra_arff, file = "oral.arff")
