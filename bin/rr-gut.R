#!/usr/bin/env Rscript

set.seed(0451)
args = commandArgs(trailingOnly=TRUE)

gut <- readRDS(args[[1]])
gut <- phyloseq::transform_sample_counts(gut, function(x) x / sum(x))
phyloseq::sample_data(gut)$cohort <-
  factor(phyloseq::sample_data(gut)$cohort,
         levels = c("Control", "Remission", "MDD"))
gut_y <- as.numeric(phyloseq::sample_data(gut)$cohort)

# csv file ---------------------------------------------------------------------
# objects: rows, attributes: columns
ra_annot <-
  cbind(data.frame(
    phyloseq::otu_table(gut),
    Class = gut_y,
    check.names = FALSE
  ))
write.table(
  ra_annot,
  file = "gut.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# arff file --------------------------------------------------------------------

otu <- data.frame(phyloseq::otu_table(gut))

# collapse a phyloseq taxonomy table to a semicolon separated string
parse_taxonomy <- function(x) {
  tax <-
    data.frame(
      tax = apply(phyloseq::tax_table(x), 1, paste, collapse = ";"),
      stringsAsFactors = FALSE
    )
  return(tibble::rownames_to_column(tax, "seq"))
}

tax <- parse_taxonomy(gut)
tax$tax <- make.unique(tax$tax, sep = "_")
write.table(tax, file = "taxonomy_key.txt", sep = "\t", quote = FALSE)

# replace colnames (sequence variants) with a taxonomy
# makes rules easier to interpret 
colnames(otu) <- tax$tax[match(colnames(otu), tax$seq)]

ra_arff <-
  cbind(
    data.frame(
      otu,
      Class = phyloseq::sample_data(gut)$cohort,
      check.names = FALSE
    )
  )

foreign::write.arff(ra_arff, file = "gut.arff")
