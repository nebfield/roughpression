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

otu <- data.frame(phyloseq::otu_table(oral))

# collapse a phyloseq taxonomy table to a semicolon separated string
parse_taxonomy <- function(x) {
  tax <-
    data.frame(
      tax = apply(phyloseq::tax_table(x), 1, paste, collapse = ";"),
      stringsAsFactors = FALSE
    )
  return(tibble::rownames_to_column(tax, "seq"))
}

tax <- parse_taxonomy(oral)
tax$tax <- make.unique(tax$tax, sep = "_")
write.table(tax, file = "taxonomy_key.txt", sep = "\t", quote = FALSE)

# replace colnames (sequence variants) with a taxonomy
# makes rules easier to interpret 
colnames(otu) <- tax$tax[match(colnames(otu), tax$seq)]

ra_arff <-
  cbind(
    data.frame(
      otu,
      Class = phyloseq::sample_data(oral)$cohort,
      check.names = FALSE
    )
  )

foreign::write.arff(ra_arff, file = "oral.arff")
