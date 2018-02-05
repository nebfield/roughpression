#!/usr/bin/env Rscript

set.seed(0451)
args = commandArgs(trailingOnly=TRUE)

gut <- readRDS(args[[1]])
gut <- phyloseq::transform_sample_counts(gut, function(x) x / sum(x))
phyloseq::sample_data(gut)$cohort <-
  factor(phyloseq::sample_data(gut)$cohort,
         levels = c("Control", "Remission", "MDD"))
gut_y <- as.numeric(phyloseq::sample_data(gut)$cohort)

# objects: rows, attributes: columns
ra_annot <- cbind(data.frame(phyloseq::otu_table(gut), Class = gut_y, check.names = FALSE))
foreign::write.arff(ra_annot, file = "gut.arff")
write.table(ra_annot, file = "gut.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)

# java -jar mahout-extensions-standalone-reducts.jar -i myData.csv -numSub 5000 -subCard 150
# java -jar mahout-extensions-standalone-reducts.jar -i gut.csv -numSub 5000 -subCard 10 -seed 0451 >test2.txt
