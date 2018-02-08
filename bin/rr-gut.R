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
ra_arff <-
  cbind(
    data.frame(
      phyloseq::otu_table(gut),
      Class = phyloseq::sample_data(gut)$cohort,
      check.names = FALSE
    )
  )

short_names <-
  tibble::tibble(short = make.unique(stringr::str_sub(
    colnames(ra_arff), start = 1, end = 5)), 
  long = colnames(ra_arff))

colnames(ra_arff) <- short_names$shorts

write.table(
  short_names,
  file = "short_names.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

foreign::write.arff(ra_arff, file = "gut.arff")

# java -jar mahout-extensions-standalone-reducts.jar -i myData.csv -numSub 5000 -subCard 150
# java -jar mahout-extensions-standalone-reducts.jar -i gut.csv -numSub 5000 -subCard 10 -seed 0451 >test2.txt